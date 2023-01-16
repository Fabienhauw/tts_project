clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC (copy)', 'FLC', 'FLC_func'});
S(mask) = [];
nb_subj=numel(S);
subjs={};
path_to_scan = {};
path_to_all_stats = {};
counter = 0;
clear a
dcm_flc_cr_matrix_design;
temp_a=a;
clear a

for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/rosso_controls';
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

for k = 1:35
    path_to_stats = path_to_all_stats{k};
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    cd (res_path)
    all_DCM = dir('DCM*');
    path_to_dcm = {};    
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    clear a
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL COMPARISON           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %% To compare models with more modulation constraints; modulations allowed on 2 connexions (stgm:mtg/smg); This is the important part for FLC
    
    clear matlabbatch
    if ~isdir(fullfile(res_path,'less_modul'))
        mkdir(fullfile(res_path,'less_modul'));
    end
    matlabbatch{1}.spm.dcm.bms.inference.dir = {fullfile(res_path,'less_modul')};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_with_new_modul;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    
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
    for fam=1:numel(path_to_dcm_with_new_modul)
        expr        = 'struct_\d';
        structur    = regexp(path_to_dcm_with_new_modul(fam), expr,'match');
        expr        = '\d';
        families    = regexp(structur{1}, expr,'match');
        fam_num     = families{1};
        partition(fam)   = str2double(fam_num{1});
    end
    family.partition = partition;
    save(fullfile(res_path, 'family_modul_constraints.mat'), 'family');
    
    if ~isdir(fullfile(res_path,'less_modul/fam_comparison'))
        mkdir(fullfile(res_path,'less_modul/fam_comparison'))
    end
    
    matlabbatch{1}.spm.dcm.bms.inference.dir = {fullfile(res_path,'less_modul/fam_comparison')};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_with_new_modul;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {fullfile(res_path,'family_modul_constraints.mat')};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    spm_jobman('run', matlabbatch)
    
    % then identify the winning family:
    BMS = load(fullfile(res_path,'less_modul/fam_comparison/BMS.mat'));
    BMS = BMS.BMS;
    [val, idx]      = max(BMS.DCM.ffx.family.post);
    win_fam         = sprintf('struct_%d',idx);
    res_path_fam    = fullfile(res_path,'less_modul',win_fam);
    if ~isdir(res_path_fam)
        mkdir(res_path_fam)
    end
    
    %identify the DCM within the winning family:
    DCM_winning_fam     = sprintf('%s',win_fam);
    common_DCM          = ~cellfun(@isempty,regexp(path_to_dcm_with_new_modul,DCM_winning_fam));
    DCM_winning_fam     = path_to_dcm_with_new_modul(common_DCM);
    path_to_dcm_lessmodul_win_fam = {};
    for dcm=1:length(DCM_winning_fam)
        path_to_dcm_lessmodul_win_fam= [path_to_dcm_lessmodul_win_fam;DCM_winning_fam(dcm)];
    end
    
    %new DCM comparison:
    clear matlabbatch
    matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_lessmodul_win_fam;
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
    %identify the DCM within the winning family:
    DCM_fam_9 = sprintf('struct_9');
    common_DCM          = ~cellfun(@isempty,regexp(path_to_dcm_with_new_modul,DCM_fam_9));
    DCM_fam_9     = path_to_dcm_with_new_modul(common_DCM);
    
    for dcm=1:length(DCM_fam_9)
        path_to_dcm_fam_9 = [path_to_dcm_fam_9;DCM_fam_9(dcm)];
    end
    
    res_path_fam_9 = fullfile(res_path, 'less_modul/struc_9');
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
