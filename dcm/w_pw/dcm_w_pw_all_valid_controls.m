% This batch script is to create a DCM model for flc for the run words vs
% pw
clear;clc;
addpath('/home/fabien.hauw/MATLAB Add-Ons');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers/DCM/w_pw');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   uncomment to apply the scripts to all subjects/controls   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..', 'FLC_func'});
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

% to know which VOI are empty, and exclude subjects with empty VOIs.
for k = 1:nb_subj
    if exist(fullfile(path_to_scan{k}, 'Aud/stats/run/stc'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/stc/dcm_model_param_modul_endroit');
    else
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    end
%     path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model');
    cd(path_to_stats)
    voi=dir('VOI_dcm*.mat');
    for v=1:length(voi)
        load(voi(v).name);
%         voi_dcm{k+1,v+1}=xY.xyz;
        if k==1
            voi_dcm{1,v+1}=xY.name;
            voi_dcm{k+1,v+1}=xY.xyz;
        else
            for dcm=1:size(voi_dcm,2)
                if strcmp(xY.name,voi_dcm{1,dcm})
                    voi_dcm{k+1,dcm}=xY.xyz;
                end
            end
        end
        if v==1
            voi_dcm{k+1,1}=subjs(k);
        end
    end
end;
for k=1:nb_subj
    if sum(cellfun(@isempty,voi_dcm(k+1,2:end)))>0
        mask(k)=1;
    else
        mask(k)=0;
    end
end

new_nsubj = nb_subj - sum(mask);
path_to_scan(mask) = '';
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
dcm_flc_w_pw_matrix_design;
temp_a=a;
clear a;

% Initialise SPM
%--------------------------------------------------------------------------
% spm('Defaults','fMRI');
% spm_jobman('initcfg');
%spm_get_defaults('cmdline',1);

for k=2:new_nsubj
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    cd (fullfile(path_to_scan{k}, 'Aud/nb_vol/run'));
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'EchoTime', 'RepetitionTime'}, {'num', 'num'});
    TE          = res{1}/1000; % ms -> seconds
    TR          = res{2}/1000; % ms -> second
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       DYNAMIC CAUSAL MODELLING                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear DCM
    
    load(fullfile(path_to_stats,'SPM.mat'));
    
    reg = {'stgm', 'smg', 'mtg mvpa', 'vwfa'};
    nregions    = length(reg);
    
    stgm=1; smg=2; mtg=3; vwfa=4;
    
    f = {
        fullfile(path_to_stats,'VOI_dcm_STGm_perso_1.mat');
        fullfile(path_to_stats,'VOI_dcm_SMG_1.mat');
        fullfile(path_to_stats,'VOI_dcm_mtg_mvpa_1.mat');
        fullfile(path_to_stats,'VOI_dcm_vwfa_sphere_1.mat')
        };
    
    for r = 1:length(f)
        XY = load(f{r});
        xY(r) = XY.xY;
    end
    
    res_path = fullfile(path_to_stats, 'dcm/all_models');
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
%     prev_dcm = dir ('DCM*');
%     if ~isempty(prev_dcm)
%         delete(prev_dcm(:).name)
%     end

    
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
    
    a = temp_a{1,9};
    b(:,:,words) = temp_a{57,9};
    mod_name = sprintf('all_models_struct_9_num_56');
    
    s.a=a;
    s.b=b;
    s.name=mod_name;
    DCM = spm_dcm_specify(SPM,xY,s);
    
    cd(res_path)
    all_DCM = dir('DCM*');
    path_to_dcm = {};
    
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    clear a
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL ESTIMATION           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %DCM Estimation
    %     clear matlabbatch
    matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = path_to_dcm;
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    % to look at results: spm_dcm_fmri_check
end