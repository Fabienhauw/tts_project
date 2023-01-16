clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/w_pw');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12')) 

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'FLC (copy)', 'FLC', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
counter = 0;
clear a
dcm_flc_w_pw_matrix_design;
temp_a=a;
clear a

for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end

for k = 1:22
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    
    cd(res_path);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL COMPARISON           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    all_DCM = dir('DCM*');
    path_to_dcm = {};
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    %to compare all DCM models;
    clear matlabbatch
    matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
    
    % fam comparison
    family.names={};
    names_fam = '';
    partition = [];
    families    = [];
    family.partition='';
    for fam=1:size(temp_a,2)
        names_fam{fam} = sprintf('F%d',fam);
    end
    family.names=names_fam;
    for fam=1:numel(path_to_dcm)
        expr        = 'struct_\d';
        structur    = regexp(path_to_dcm(fam), expr,'match');
        expr        = '\d';
        families    = regexp(structur{1}, expr,'match');
        fam_num     = families{1};
        partition(fam)   = str2double(fam_num{1});
    end
    family.partition = partition;
    save('family.mat', 'family');
    
    if ~isdir(fullfile(res_path,'fam_comparison'))
        mkdir(fullfile(res_path,'fam_comparison'))
    end
    
    matlabbatch{1}.spm.dcm.bms.inference.dir = {fullfile(res_path,'fam_comparison')};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {fullfile(res_path,'family.mat')};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
    % then identify the winning family:
    BMS = load(fullfile(res_path,'fam_comparison/BMS.mat'));
    BMS = BMS.BMS;
    [val, idx] = max(BMS.DCM.ffx.family.post);
    win_fam         = sprintf('struct_%d',idx);
    res_path_fam = fullfile(res_path, win_fam);
    if ~isdir(res_path_fam)
        mkdir(res_path_fam)
    end
    
    %identify the DCM within the winning family:
    DCM_winning_fam = dir(sprintf('*%s*.mat',win_fam));
    path_to_dcm_win_fam = {};
    for dcm=1:length(DCM_winning_fam)
        new_dcm      = fullfile(res_path,DCM_winning_fam(dcm).name);
        path_to_dcm_win_fam= [path_to_dcm_win_fam;new_dcm];
    end
    
    %new DCM comparison:
    matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_win_fam;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    
    
    %% to compare systematically all models within the best group model (ie
    % struct 9)
    
    path_to_dcm_fam_9 = {};
    %identify the DCM within the 9th family:
    DCM_fam_9 = dir(sprintf('*struct_9*.mat'));
    for dcm=1:length(DCM_fam_9)
        new_dcm             = fullfile(res_path,DCM_fam_9(dcm).name);
        path_to_dcm_fam_9   = [path_to_dcm_fam_9;new_dcm];
    end
    
    res_path_fam_9 = fullfile(res_path, 'struc_9');
    if ~isdir(res_path_fam_9)
        mkdir(res_path_fam_9)
    end

    %new DCM comparison:
    matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam_9};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_fam_9;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    clear matlabbatch
end