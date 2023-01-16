clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/w_pw');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12')) 

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'FLC', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
counter = 0;

for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end
for k = 1:22 %parfor 1:numel(S)
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path)
    all_DCM = dir('DCM*');
    all_DCM_with_modul = all_DCM;
    path_to_dcm = {};
    path_to_dcm_with_modul = {};
    DCM_without_modul = dir('DCM*num_1_*');
    mask = ismember({all_DCM_with_modul.name}, {DCM_without_modul.name});
    all_DCM_with_modul(mask) = [];
    for dcm = 1 : length(all_DCM_with_modul)
        new_dcm      			= fullfile(res_path,all_DCM_with_modul(dcm).name);
        path_to_dcm_with_modul  	= [path_to_dcm_with_modul;new_dcm];
    end
    
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL ESTIMATION           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %DCM Estimation
    %     clear matlabbatch
    for dcm_mod = 1:length(path_to_dcm)
        counter = counter + 1;
        
%         matlabbatch{counter}.spm.dcm.estimate.dcms.subj(1).dcmmat = {''};
%         matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
%         matlabbatch{counter}.spm.dcm.estimate.est_type = 3;
%         matlabbatch{counter}.spm.dcm.estimate.fmri.analysis = 'time';

        matlabbatch{counter}.spm.dcm.estimate.dcms.subj.dcmmat = path_to_dcm(dcm_mod);
        matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
    end
end

cd('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/w_pw')

par.run = 1;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_wpw_estimate';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);
