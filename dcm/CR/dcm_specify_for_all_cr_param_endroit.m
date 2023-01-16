% This batch script is to create a DCM model for flc for the CR
% better use this script than the "reduced one", because this simulates all
% models (longer) with all possible modulations, but in the end only takes into
% account the ones with "new reduced" modulation...

clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12')) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC', 'FLC_func'});
S(mask) = [];
% mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC_func'});

subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats = {};
repet_time = {};
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
    repet_time          = [repet_time;1.3];
end

D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/rosso_controls';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'README'});
S(mask) = [];

nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    repet_time          = [repet_time;3];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
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

for k = 18 % parfor k = 1 : numel(S)
    matlabbatch=[];
    b=[];
    family=[];
    path_to_stats = path_to_all_stats{k};
    
    cd (fullfile(path_to_scan{k}, 'Aud/nb_vol/chaperon'));
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'EchoTime', 'RepetitionTime'}, {'num', 'num'});
    TE          = res{1}/1000; % ms -> seconds
    TR          = res{2}/1000; % ms -> second
    
    
    % Initialise SPM
    %--------------------------------------------------------------------------
    %spm('Defaults','fMRI');
    %spm_jobman('initcfg');
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
        fullfile(path_to_stats,'VOI_dcm_STGm_perso_1.mat');
        fullfile(path_to_stats,'VOI_dcm_SMG_1.mat');
        fullfile(path_to_stats,'VOI_dcm_mtg_mvpa_1.mat');
        fullfile(path_to_stats,'VOI_dcm_vwfa_sphere_1.mat')
        };
    
    for r = 1:4
       XY = load(f{r});
       xY(r) = XY.xY;
    end
    
    res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models');
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    dcm_flc_cr_matrix_design;
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
    
    name_count=0;
    length(find(~cellfun(@isempty,temp_a(2:end,:))));
    for i=1:size(temp_a,2)
        a=temp_a{1,i};
        modul_models = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j=2:modul_models
            if i>1
                name_count = length(find(~cellfun(@isempty,temp_a(2:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(2:j,i))));
            else
                name_count = length(find(~cellfun(@isempty,temp_a(2:j,i))));
            end
            b(:,:,words)=temp_a{j,i};
            mod_name = sprintf('all_models_struct_%d_num_%d_tot_%d',i,j-1,name_count);
            s.a=a;
            s.b=b;
            s.name=mod_name;
            %%%%%%%%this line below to comment to avoid re defining
            %%%%%%%%models
            DCM = spm_dcm_specify(SPM,xY,s);
        end
    end
end