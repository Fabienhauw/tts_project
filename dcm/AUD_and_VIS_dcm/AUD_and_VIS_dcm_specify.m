% This batch script is to create a DCM model for every participant

clear;clc;
addpath(genpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm'));
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sphere_radius = 4;
model_kind = 3;
if model_kind == 1
    dcm_folder = 'dcm_model_param_modul_speech_baseline';
elseif model_kind == 2
    dcm_folder = 'dcm_model_param_modul_sent_scramble';
    elseif model_kind == 3
    dcm_folder = 'dcm_model_param_modul';
end

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_syn = S;
S_syn(mask_gauch) = [];

S = [S_syn ; S_con];

subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats_aud = {};
path_to_all_stats_vis = {};
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path_aud            = fullfile(path_to_subj, 'Aud/loc/stats_s5_without_resting');
    tmp_path_vis            = fullfile(path_to_subj, 'Vis/loc/stats_s5_without_resting');
    path_to_all_stats_aud   = [path_to_all_stats_aud; tmp_path_aud];
    path_to_all_stats_vis   = [path_to_all_stats_vis; tmp_path_vis];
end

%% model parameters
sounds = 1; speech = 2;

subj={subjs};
nsubjects   = length(subj);

reg         = {'lSTG_common', 'SMG_common', 'VWFA_common'};
ROIs_coord = [-58 -11 0; -48 -41 16; -42 -38 -20];

nregions    = length(reg);
lpstg = 1; smg = 2; vwfa = 3;

% this code is to check how much of the variance in your VOI is explained by the timeseries :
% load('VOI_****_adapted_1.mat')
% 100*xY.s(1)/sum(xY.s)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory containing the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AUD_and_VIS_dcm_matrix_3roi_design;

temp_a=a;
clear a

for k = 1 : numel(S)
    matlabbatch=[];
    b=[];
    family=[];
    path_to_stats = path_to_all_stats_aud{k};
    
    cd (fullfile(path_to_scan{k}, 'Aud/loc/param'));
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'EchoTime', 'RepetitionTime'}, {'vect', 'num'});
    res{1} = unique(res{1}); res{1} = sort(res{1}, 'ascend');
    TE          = res{1}(2)/1000; % ms -> seconds
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
    
    SPM = load(fullfile(path_to_stats,dcm_folder,'SPM.mat'));
    SPM = SPM.SPM;
    
    reg = {'lSTG_common', 'SMG_common', 'VWFA_common'};
    nregions    = length(reg);
    
    lpstg = 1; smg = 2; vwfa = 3;
    
    f = {
        fullfile(path_to_stats,sprintf('VOI_lSTG_common_-58_-11_0_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_SMG_common_-48_-41_16_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_VWFA_common_-42_-38_-20_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        };
    
    cd (path_to_stats)
    for r = 1 : nregions
       XY = load(f{r});
       xY(r) = XY.xY;
    end
    
    res_path = fullfile(path_to_stats, dcm_folder, sprintf('all_3_rois_models_%dmm', sphere_radius));
    
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    
    %% common parameters for all models;
    include = [1 1 0];
    
    % C-matrix = driving / input matrice (lines for roi, columns for conditions)
    c = zeros(nregions, sum(include)); % better use sum of include than nconditions, because if you don't include some conditions, it will bug the estimation step;
    c(lpstg,sounds)   = 1;
    % D-matrix (disabled)
    d = zeros(nregions,nregions,0);

    
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
        a=temp_a{1,i}(:,:,1);
        modul_models = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j=1:modul_models
            if i>1
                name_count = length(find(~cellfun(@isempty,temp_a(1:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(1:j,i))));
            else
                name_count = length(find(~cellfun(@isempty,temp_a(1:j,i))));
            end
            
            b(:,:,speech)=temp_a{j,i}(:,:,2);
            
            mod_name = sprintf('struct_%d_num_%d_tot_%d',i,j,name_count);
            s.a=a;
            s.b=b;
            s.name=mod_name;
            %%%%%%%%this line below to comment to avoid re defining
            %%%%%%%%models
            DCM = spm_dcm_specify(SPM,xY,s);
        end
    end
end

%% same for visual models
clear a
AUD_and_VIS_dcm_matrix_3roi_design;

temp_a=a;
clear a

for k = 1 : numel(S)
    matlabbatch=[];
    b=[];
    family=[];
    path_to_stats = path_to_all_stats_vis{k};
    
    cd (fullfile(path_to_scan{k}, 'Vis/loc/param'));
    json=dir('*.json');
    json=json.name;
    
    res = get_string_from_json(json, {'EchoTime', 'RepetitionTime'}, {'num', 'num'});
    res{1} = unique(res{1}); res{1} = sort(res{1}, 'ascend');
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
    
    SPM = load(fullfile(path_to_stats,'dcm_model_param_modul','SPM.mat'));
    SPM = SPM.SPM;
    
    reg = {'common_occip', 'common_VWFA', 'common_SMG'};
    nregions    = length(reg);
    
    lpstg = 1; smg = 2; vwfa = 3;
    
    f = {
        fullfile(path_to_stats,sprintf('VOI_common_occip_-28_-84_-12_%dmm_sph_DCM_ROI_adapted_to_vis_con12_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_common_VWFA_-50_-54_-14_%dmm_sph_DCM_ROI_adapted_to_vis_con12_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_common_SMG_-52_-44_20_%dmm_sph_DCM_ROI_adapted_to_vis_con12_adj_eoi5_1.mat', sphere_radius));
        };
    
    cd (path_to_stats)
    for r = 1 : nregions
       XY = load(f{r});
       xY(r) = XY.xY;
    end
    
    res_path = fullfile(path_to_stats, 'dcm_model_param_modul', sprintf('all_3_rois_models_%dmm', sphere_radius));
    
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    
    %% common parameters for all models;
    include = [1 1];
    
    % C-matrix = driving / input matrice (lines for roi, columns for conditions)
    c = zeros(nregions, sum(include)); % better use sum of include than nconditions, because if you don't include some conditions, it will bug the estimation step;
    c(lpstg,sounds)   = 1;
    % D-matrix (disabled)
    d = zeros(nregions,nregions,0);

    
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
        a=temp_a{1,i}(:,:,1);
        modul_models = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j=1:modul_models
            if i>1
                name_count = length(find(~cellfun(@isempty,temp_a(1:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(1:j,i))));
            else
                name_count = length(find(~cellfun(@isempty,temp_a(1:j,i))));
            end
            
            b(:,:,speech)=temp_a{j,i}(:,:,2);
            
            mod_name = sprintf('struct_%d_num_%d_tot_%d',i,j,name_count);
            s.a=a;
            s.b=b;
            s.name=mod_name;
            %%%%%%%%this line below to comment to avoid re defining
            %%%%%%%%models
            DCM = spm_dcm_specify(SPM,xY,s);
        end
    end
end