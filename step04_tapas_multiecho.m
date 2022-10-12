% #TAPAS/PhysIO -> SPM12 (https://github.com/translationalneuromodeling/tapas)

clear
clc

load e

addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/abin'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB'));

CLUSTER = 0; % autoaddobjet = 0 when par.sge = 1; so no cluster for this one;


%% fetch physio files


main_dir_physio = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';
subj_to_preproc = 'Control18';
% subj_to_preproc = 'Sujet\d[4-9]|Sujet1[0-2]|Sujet1[4-9]|Sujet[2][0-9]';
% subj_to_preproc = 'Control\d[2-9]|Control[1][0-2]';

e_physio = exam(main_dir_physio, subj_to_preproc); % all subjects with multi-echo

% beware because order of the sequences matters here, while in step01, it is
% automatically reordered according to the num of the sequence (S03, S04,
% etc)
e_physio.addSerie('RS_PhysioLog$','physio_rs',1)
e_physio.addSerie('visual_PhysioLog$','physio_visual',1)
e_physio.addSerie('unframed_color_PhysioLog$','physio_unframed_color',1)
% % e_physio.addSerie('color_PhysioLog$','physio_color',1)
e_physio.addSerie( 'audio_PhysioLog$','physio_audio' ,1)


e_physio.getSerie('physio').addVolume('dcm$','dcm',1)

e_physio.getSerie('physio').addPhysio('Info.log$','info',1)
e_physio.getSerie('physio').addPhysio('PULS.log$','puls',1)
e_physio.getSerie('physio').addPhysio('RESP.log$','resp',1)

% e_physio.explore

% e_physio.getSerie('physio').getVolume('info').removeEmpty
% e_physio.getSerie('physio').getVolume('puls').removeEmpty
% e_physio.getSerie('physio').getVolume('resp').removeEmpty

info = e_physio.getSerie('physio').getPhysio('info').removeEmpty;
puls = e_physio.getSerie('physio').getPhysio('puls').removeEmpty;
resp = e_physio.getSerie('physio').getPhysio('resp').removeEmpty;

% e.getSerie('run').getRP('spm').plot()
%% Prepare dirs & files
run = e.getSerie('run');

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
job_afni_remove_nan( run.getVolume('^wts_OC'), par );


%%

volume = run.getVolume('^nwts_OC').removeEmpty;

outdir = get_parent_path(volume.getPath);

run.addRP('tedana','rp_spm.txt','rp_spm',1)
rp     = run.getRP('rp_spm').removeEmpty;

mask   = volume.getExam.getSerie('anat').getVolume('^rwp[23]');
% mask   = volume.getExam.getSerie('anat').getVolume('^rwp[2]');


%%

clear par

%----------------------------------------------------------------------------------------------------------------------------------------------------
% ALWAYS MANDATORY
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio   = 1; % put to 0 if error during recording of physio
par.noiseROI = 1;
par.rp       = 1;

par.TR     = 1.660;
par.nSlice = 60;

par.volume = volume;
par.outdir = outdir;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Physio
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio_Info = info;
par.physio_PULS = puls;
par.physio_PULSE = puls;
par.physio_RESP = resp;

par.physio_RETROICOR        = 1;
par.physio_HRV              = 1;
par.physio_RVT              = 1;
par.logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
par.physio_logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
par.physio_slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.


%----------------------------------------------------------------------------------------------------------------------------------------------------
% noiseROI
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.noiseROI_mask   = mask;
par.noiseROI_volume = volume;

par.noiseROI_thresholds   = [0.95 0.7];     % keep voxels with tissu probabilty >= 95%
par.noiseROI_n_voxel_crop = [2 1];           % crop n voxels in each direction, to avoid partial volume
par.noiseROI_n_components = 10;              % keep n PCA componenets
% 
% par.noiseROI_thresholds   = [0.95];     % keep voxels with tissu probabilty >= 95%
% par.noiseROI_n_voxel_crop = [2];           % crop n voxels in each direction, to avoid partial volume
% par.noiseROI_n_components = 10;              % keep n PCA componenets

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Realignment Parameters
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.rp_file = rp;

par.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
par.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'
par.rp_threshold = 0.5;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Other
%----------------------------------------------------------------------------------------------------------------------------------------------------
par.print_figures = 0; % 0 , 1 , 2 , 3

% classic matvol
if CLUSTER
    par.run = 0;
    par.sge = 1;
    par.sge_queu = 'normal,bigmem';
else
    par.run = 1;
    par.sge = 0;
end
par.display  = 0;
par.redo     = 0;

% cluster
par.jobname  = 'spm_physio';
par.walltime = '04:00:00';
par.mem      = '4G';

job_physio_tapas( par );

% if error : plot((pulse.UUID-pulse.UUID(1))*0.0025, pulse.x5b7bbc23_25e7_4fd4_9dc1_ffdcbefdcb34)
% do plot physio for respi and pulse and check whether it is the same last.