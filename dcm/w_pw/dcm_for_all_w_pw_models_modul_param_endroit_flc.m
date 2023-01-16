% This batch script is to create a DCM model for flc for the run words vs
% pw
% Resp path is called 'minus Sujet08, but this subject is included in interindividual DCM comparison,
% only excluded from CR DCM analysis and ROI creation.
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
repet_time = {};
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
    repet_time = [repet_time;1.3];
end

%% model parameters
nconditions = 2;

sounds=1; words=2;

subj={subjs};
nsubjects   = length(subj);

reg         = {'stgm', 'smg', 'mtg mvpa', 'vwfa'};
nregions    = length(reg);

stgm=1; smg=2; mtg=3; vwfa=4;
% this code is to check how much of the variance in your VOI is explained by the timeseries :
% load('VOI_****_adapted_1.mat')
% 100*xY.s(1)/sum(xY.s)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory containing the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 20 %parfor 1:numel(S)
    
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    new_path = fullfile(path_to_stats, 'dcm_minus_Sujet08');
    if ~isdir(new_path)
        mkdir(new_path)
    end
    cd (new_path)
    if isempty(dir('VOI_dcm_STGm_perso_1*'))
        copyfile(fullfile(path_to_stats,'VOI_dcm_STGm_perso_1.mat'), fullfile(new_path,'VOI_dcm_STGm_perso_1.mat'));
        copyfile(fullfile(path_to_stats,'VOI_dcm_SMG_1.mat'), fullfile(new_path,'VOI_dcm_SMG_1.mat'));
        copyfile(fullfile(path_to_stats,'VOI_dcm_mtg_mvpa_1.mat'), fullfile(new_path,'VOI_dcm_mtg_mvpa_1.mat'));
        copyfile(fullfile(path_to_stats,'VOI_dcm_vwfa_sphere_1.mat'), fullfile(new_path,'VOI_dcm_vwfa_sphere_1.mat'));
    end
    cd (fullfile(path_to_scan{k}, 'Aud/nb_vol/run'));
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'EchoTime', 'RepetitionTime'}, {'num', 'num'});
    TE          = res{1}/1000; % ms -> seconds
    TR          = res{2}/1000; % ms -> second
    
    
    % Initialise SPM
    %--------------------------------------------------------------------------
    % spm('Defaults','fMRI');
    % spm_jobman('initcfg');
    %spm_get_defaults('cmdline',1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       DYNAMIC CAUSAL MODELLING                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear DCM
    
    SPM = load(fullfile(path_to_stats,'SPM.mat'));
    SPM = SPM.SPM;
    
    reg = {'stgm', 'smg', 'mtg mvpa', 'vwfa'};
    nregions    = length(reg);
    
    stgm=1; smg=2; mtg=3; vwfa=4;
    
    f = {
        fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_STGm_perso_1.mat');
        fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_SMG_1.mat');
        fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_mtg_mvpa_1.mat');
        fullfile(path_to_stats,'dcm_minus_Sujet08/VOI_dcm_vwfa_sphere_1.mat')
        };
    
    clear xY
    tmp_xY = [];
    for r = 1:length(reg)
        XY = load(f{r});
        tmp_xY = [tmp_xY;XY];
    end
    xY = [tmp_xY.xY];
    
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    dcm_flc_w_pw_matrix_design;
    temp_a=a;
    clear a
    
    %% common parameters for all models;
    % C-matrix = driving
    c = zeros(nregions, nconditions);
    c(stgm,sounds)   = 1;
    % D-matrix (disabled)
    d = zeros(nregions,nregions,0);
    
    include = [1 1];
    
    s = struct();
    s.name       = 'to change';
    s.u          = include;                 % Conditions
    s.delays     = repmat(TR/2,1,nregions)';   % Slice timing for each region
    s.TE         = TE;
    s.nonlinear  = false;
    s.two_state  = false;
    s.stochastic = false;
    s.centre     = false;
    s.induced    = 0;
    s.c          = c;
    s.d          = d;
    
    name_count = 0;
    length(find(~cellfun(@isempty,temp_a(2:end,:))))
    for i = 1:size(temp_a,2)
        a                   = temp_a{1,i};
        modul_models        = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j = 2:modul_models
            if i>1
                name_count  = length(find(~cellfun(@isempty,temp_a(2:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(2:j,i))));
            else
                name_count  = length(find(~cellfun(@isempty,temp_a(2:j,i))));
            end
            b(:,:,words)    = temp_a{j,i};
            mod_name        = sprintf('all_models_struct_%d_num_%d_tot_%d',i,j-1,name_count);
            s.a             = a;
            ss.b            = b;
            s.name          = mod_name;
            
            %%%%%%%%this line below to comment to avoid re defining
            %%%%%%%%models
%             DCM = spm_dcm_specify(SPM,xY,s);
        end
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
    matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = path_to_dcm;
    
    
    %%%%%%%% this line below to comment to avoid re estimating
    %%%%%%%% models
    
%     spm_jobman('run',matlabbatch);
    clear matlabbatch
    % to look at results: spm_dcm_fmri_check
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL COMPARISON           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    clear matlabbatch
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
end

% %to compare only DCM models with modulation;
% clear matlabbatch
% res_path_2 = fullfile(res_path,'comparison_DCM_with_modul');
% if ~isdir(res_path_2)
%     mkdir(res_path_2)
% end
% matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_2};
% matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_with_modul;
% matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
% matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
% matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
% matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
% matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
% matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
%
% spm_jobman('run', matlabbatch)
% clear matlabbatch

% max(BMS.DCM.ffx.model.post)
% find(BMS.DCM.ffx.model.post>0.3)
% find(BMS.DCM.ffx.model.post>0.01)clear;clc;
