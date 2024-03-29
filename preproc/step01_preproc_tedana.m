%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   MULTI - ECHO    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))
addpath(genpath('/network/lustre/iss02/apps/software/scit/AFNI'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/abin'))

main_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';

% subj_to_preproc = 'Sujet\d[4-9]|Sujet[1-2][0-9]|Control\d[2-9]|Control[1-2][0-2]'; % this one is for subjects who had the unframed color sequence
% subj_to_preproc = 'Sujet0[123]|Control01'; % this one is for subjects who had the framed color sequence;
subj_to_preproc = 'Sujet01';
% subj_to_preproc = '';
e = exam(main_dir, subj_to_preproc); % all subjects with multi-echo

% e = exam(main_dir, 'SYNESTHEX_Sujet'); % all subjects with multi-echo

%% Get files paths #matvol

% % Anat
e.addSerie('T1w$', 'anat_T1', 1 );
e.getSerie('anat').addVolume('^v_.*nii','v',1);

% Func
% e.addSerie('Localizer_visual$' , 'run_localizer_visual', 1 );
% e.addSerie('Localizer_audio$'  , 'run_localizer_audio' , 1 );
e.addSerie('Localizer_unframed_color$'  , 'run_localizer_unframed_color' , 1 );
% e.addSerie('Localizer_color$'  , 'run_localizer_color' , 1 );
% e.addSerie('RS$'  , 'run_rs' , 1 );

e.getSerie('run').addVolume('^v_.*nii','v',3);

e.reorderSeries('name')

e.explore


%% Cluster ?

CLUSTER = 1;


%% segment cat12 #CAT12->SPM12

anat = e.gser('anat_T1').gvol('^v');

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
    par.mem      = '8G'; 
else
    par.run = 1;
    par.sge = 0;
end
par.redo    = 0;
par.display = 0;

par.subfolder = 0;

par.GM        = [1 1 1 0]; % warped_space_Unmodulated (wp1*)     / warped_space_modulated (mwp1*)     / native_space (p1*)     / native_space_dartel_import (rp1*)
par.WM        = [1 1 1 0]; %                          (wp2*)     /                        (mwp2*)     /              (p2*)     /                            (rp2*)
par.CSF       = [1 1 1 0]; %                          (wp3*)     /                        (mwp3*)     /              (p3*)     /                            (rp3*)
par.TPMC      = [0 0 0 0]; %                          (wp[456]*) /                        (mwp[456]*) /              (p[456]*) /                            (rp[456]*)   This will create other probalities map (p4 p5 p6)

par.label     = [1 1 0] ;  % native (p0*)  / normalize (wp0*)  / dartel (rp0*)       This will create a label map : p0 = (1 x p1) + (3 x p2) + (1 x p3)
par.bias      = [1 1 0] ;  % native (ms*)  / normalize (wms*)  / dartel (rms*)       This will save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % native (mis*) / normalize (wmis*) / dartel (rmis*)       This will save the bias field corrected  + SANLM (local) T1

par.warp      = [1 1];     % warp fields  : native->template (y_*) / native<-template (iy_*)

par.doSurface = 0;
par.jacobian  = 0;         % write jacobian determinant in normalize space
par.doROI     = 0;         % will compute the volume in each atlas region

job_do_segmentCAT12(anat,par);




%% Sort echos #matvol

clear par
par.run  = 1;
par.fake = 0;
par.sge  = 0;

par.redo = 1;
meinfo = job_sort_echos( e.getSerie('run') , par );


%% job_afni_proc_multi_echo #ANFI

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake = 0;
par.redo = 1;

par.seperate = 1;
par.write_nifti = 1;

par.subdir = 'afni_vtd';
par.blocks  = {'despike', 'tshift', 'volreg'};
job_afni_proc_multi_echo( meinfo, par );


%% do_fsl_robust_mask_epi #FSL

fin  = e.getSerie('run').getVolume('^vtde1');

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 1;
par.fsl_output_format = 'NIFTI_GZ';
do_fsl_robust_mask_epi( fin, par );

% Checkpoint & unzip
par.jobname = 'unzip_and_keep__bet';
e.getSerie('run').getVolume('^bet_Tmean_vtde1$').removeEmpty.unzip_and_keep(par)



%% TEDANA #Python

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.fake  = 0;
par.redo  = 0;
par.pct = 0;

% cluster
par.walltime = '12:00:00';      % HH:MM:SS
par.mem      = '16G';           % ICA is very memory consuming
par.sge_nb_coeur = 2;           % I dont't know why, but 2 CPU increase the "stability" of the job on the cluster

job_tedana_009a1( meinfo, 'vtd', 'tedana009a1_vtd', 'bet_Tmean_vtde1_mask.nii.gz ', par );
% job_t2smap( meinfo, 'vtd', 't2smap_vtd', 'bet_Tmean_vtde1_mask.nii.gz ', par );

% Checkpoint & unzip
par.jobname = 'unzip_and_keep__tedana';
% par.jobname = 'unzip_and_keep__t2smap';
% e.getSerie('run').getVolume('dn_ts_OC').removeEmpty.unzip_and_keep(par)
e.getSerie('run').getVolume('ts_OC').removeEmpty.unzip_and_keep(par)


%% Coregister TEDANA outputs to Anat #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo  = 0;
par.type  = 'estimate';

src = e.getSerie('run').removeEmpty.getVolume('^bet_Tmean_vtde1$');
oth = e.getSerie('run').removeEmpty.getVolume('(^ts_OC)|(^dn_ts_OC)');
ref = e.getSerie('run').removeEmpty.getExam.getSerie('anat_T1').getVolume('^p0');

par.jobname = 'spm_coreg_epi2anat';
job_coregister(src,ref,oth,par);


%% Normalize TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 0;
par.vox = [2.5 2.5 2.5]; % IMPORTANT keep original EPI voxel size
img = e.getSerie('run').getVolume('(^ts_OC)|(^dn_ts_OC)').removeEmpty;
for i =1 :length(img)
    if ~isempty(regexp(img(i).path,'gz'))
        unix(sprintf('gunzip -f %s',  img(i).path));
    end
end
img = e.getSerie('run').getVolume('(^ts_OC)|(^dn_ts_OC)').removeEmpty;
y   = img.getExam.getSerie('anat_T1').removeEmpty.getVolume('^y');
par.jobname = 'spm_normalize_epi';
job_apply_normalize(y,img,par);

% Nomalize Tmean, used later for PhysIO Noise ROI
img = e.getSerie('run').removeEmpty.getVolume('^bet_Tmean_vtde1$');
y   = img.getExam.getSerie('anat_T1').getVolume('^y').removeEmpty;
par.jobname = 'spm_normalize_meanepi';
job_apply_normalize(y,img,par);


%%  Smooth TEDANA outputs #SPM12

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 1;

img = e.getSerie('run').getVolume('^wt.*_OC').removeEmpty;

% par.smooth   = [3 3 3];
% par.prefix   = 's3';
% job_smooth(img,par);

par.smooth   = [4 4 4];
par.prefix   = 's4';
job_smooth(img,par);

par.smooth   = [5 5 5];
par.prefix   = 's5';
job_smooth(img,par);

par.smooth   = [6 6 6];
par.prefix   = 's6';
job_smooth(img,par);

% par.smooth   = [8 8 8];
% par.prefix   = 's8';
% job_smooth(img,par);

%% coregister WM & CSF on functionnal (using the warped mean) #SPM12
% This will be used for TAPAS:PhysIO

clear par
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.redo = 1;

ref = e.getSerie('run');
ref = ref(:,1).getVolume('wbet_Tmean_vtde1');
src = e.getSerie('anat_T1').getVolume('^wp2');
oth = e.getSerie('anat_T1').getVolume('^wp3');
par.type = 'estimate_and_write';
par.jobname = 'spm_coreg_WMCSF2wEPI';
job_coregister(src,ref,oth,par);


%% Save examArray

save e e
