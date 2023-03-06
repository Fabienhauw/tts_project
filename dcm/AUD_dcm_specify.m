% This batch script is to create a DCM model for flc for the CR
% better use this script than the "reduced one", because this simulates all
% models (longer) with all possible modulations, but in the end only takes into
% account the ones with "new reduced" modulation...

clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nroi = 3;

i=0;
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
S_con_app = S;
S_con_app(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];

S = [S_droit ; S_con_app];

subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats = {};
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/loc/stats_s5_without_resting');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

%% model parameters
nconditions = 3;

sounds = 1; speech = 2; words = 3;

subj={subjs};
nsubjects   = length(subj);

if nroi == 4
    reg         = {'lstgm', 'smg', 'mfg', 'vwfa'};
    ROIs_coord = [-45 -24 8; -48 -44 23; -50 6 53; -45 -51 -10];
elseif nroi == 3
    reg         = {'lpstg', 'smg', 'vwfa'};
    ROIs_coord = [-70 -28 3;-48 -44 23; -45 -51 -10];
end

nregions    = length(reg);

if nroi == 3
    lstgm=1; smg=2; vwfa=3;
elseif nroi == 4
    lstgm=1; smg=2; mfg=3; vwfa=4;
end
% this code is to check how much of the variance in your VOI is explained by the timeseries :
% load('VOI_****_adapted_1.mat')
% 100*xY.s(1)/sum(xY.s)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory containing the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : numel(S)
    matlabbatch=[];
    b=[];
    family=[];
    path_to_stats = path_to_all_stats{k};
    
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
    
    SPM = load(fullfile(path_to_stats,'dcm_model_param_modul/SPM.mat'));
    SPM = SPM.SPM;
    
    if nroi == 3
        reg         = {'lpstg', 'smg', 'vwfa'};
    elseif nroi == 4
        reg         = {'lstgm', 'smg', 'mfg', 'vwfa'};
    end
    
    nregions    = length(reg);
    
    if nroi == 3
        lstgm=1; smg=2; vwfa=3;
    elseif nroi == 4
        lstgm=1; smg=2; mfg=3; vwfa=4;
    end
    
    if nroi == 3
        f = {
        fullfile(path_to_stats,'VOI_lSTGm_-45_-24_8_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        fullfile(path_to_stats,'VOI_SMG_-48_-44_23_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        fullfile(path_to_stats,'VOI_VWFA_-45_-51_-10_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        };
    elseif nroi == 4
        f = {
        fullfile(path_to_stats,'VOI_lSTGm_-45_-24_8_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        fullfile(path_to_stats,'VOI_SMG_-48_-44_23_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        fullfile(path_to_stats,'VOI_MFG_-50_6_53_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        fullfile(path_to_stats,'VOI_VWFA_-45_-51_-10_4mm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat');
        };
    end
    
    cd (path_to_stats)
    for r = 1 : nregions
       XY = load(f{r});
       xY(r) = XY.xY;
    end
    
    if nroi == 3
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/all_4_rois_models_4mm');
    end
    
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    if nroi == 3
        AUD_dcm_matrix_3roi_design;
    elseif nroi == 4
        AUD_dcm_matrix_4roi_design;
    end
    
    temp_a=a;
    clear a
    
    %% common parameters for all models;
    
    if nroi == 3
        include = [1 1 0];
    elseif nroi == 4
        include = [1 1 1];
    end
    
    % C-matrix = driving / input matrice (lines for roi, columns for conditions)
    c = zeros(nregions, sum(include)); % better use sum of include than nconditions, because if you don't include some conditions, it will bug the estimation step;
    c(lstgm,sounds)   = 1;
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
        a=temp_a{1,i};
        modul_models = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j=2:modul_models % skip the first line which corresponds to the fixed connections;
            if i>1
                name_count = length(find(~cellfun(@isempty,temp_a(2:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(2:j,i))));
            else
                name_count = length(find(~cellfun(@isempty,temp_a(2:j,i))));
            end
            
            b(:,:,speech)=temp_a{j,i};
            if nroi == 4
                b(:,:,words)=temp_a{j,i};
            end
            
            mod_name = sprintf('struct_%d_num_%d_tot_%d',i,j-1,name_count);
            s.a=a;
            s.b=b;
            s.name=mod_name;
            %%%%%%%%this line below to comment to avoid re defining
            %%%%%%%%models
            DCM = spm_dcm_specify(SPM,xY,s);
        end
    end
end