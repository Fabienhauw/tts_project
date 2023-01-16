clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
D = '/network/lustre/iss02/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats = {};
counter = 0;

for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

D = '/network/lustre/iss02/cohen/data/Fabien_official/flc_vs_controls/fmri/rosso_controls';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'README'});
S(mask) = [];

nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

for k = 1%:35
    path_to_stats = path_to_all_stats{k};
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    cd(res_path)
    if ~isdir(res_path)
        mkdir(res_path)
    end
    all_DCM = dir('DCM*');
    all_DCM_with_modul = all_DCM;
    path_to_dcm = {};
    path_to_dcm_with_modul = {};
    DCM_without_modul = dir('DCM*num_1_*');
    mask = ismember({all_DCM_with_modul.name}, {DCM_without_modul.name});
    all_DCM_with_modul(mask) = [];
    for dcm=1:length(all_DCM_with_modul)
        new_dcm      = fullfile(res_path,all_DCM_with_modul(dcm).name);
        path_to_dcm_with_modul  = [path_to_dcm_with_modul;new_dcm];
    end
    
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    clear a
    clear b
    
    path_to_dcm_with_new_modul = {};
    dcm_flc_cr_matrix_design_new_reduced_modul;
    for jj=1:length(all_DCM)
        DCM = load(all_DCM(jj).name);
        DCM = DCM.DCM;
        dcm_ok=DCM.b(:,:,2);
        if isempty(find(dcm_ok(B==0)==1))
            path_to_dcm_with_new_modul = [path_to_dcm_with_new_modul;fullfile(res_path,all_DCM(jj).name)];
        end
    end
    
    clear a
    clear B
    clear dcm_ok
    
    B=zeros(nbreg);
    B(smg,stgm) = 1;
    B(mtg,stgm) = 1;
    B(mtg,smg) = 1;
    B(smg,mtg) = 1;
    path_to_dcm_with_interm_modul = {};
    for jj=1:length(all_DCM)
        DCM = load(all_DCM(jj).name);
        DCM = DCM.DCM;
        dcm_ok=DCM.b(:,:,2);
        if isempty(find(dcm_ok(B==0)==1))
            path_to_dcm_with_interm_modul = [path_to_dcm_with_interm_modul;fullfile(res_path,all_DCM(jj).name)];
        end
    end
    clear a
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL ESTIMATION           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %DCM Estimation
    for dcm_mod = 1:length(path_to_dcm)
        counter = counter + 1;
        matlabbatch{counter}.spm.dcm.estimate.dcms.subj.dcmmat = path_to_dcm(dcm_mod);
        matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
    end
    %%%%%%%% this line below to comment to avoid re estimating
    %%%%%%%% models
    
end

cd('/network/lustre/iss02/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_cr_estimate';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);
