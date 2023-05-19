% This batch script is to create a DCM model for every participant

clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nroi = 4;
lexic = 1;
sphere_radius = 6;
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
sounds = 1; speech = 2; words = 3;

subj={subjs};
nsubjects   = length(subj);

if nroi == 4
    reg         = {'lstgm', 'smg', 'sts', 'vwfa'};
    ROIs_coord = [-45 -24 8; -48 -44 23; -65 -41 6; -45 -51 -10];
elseif nroi == 3
    reg         = {'lstgm', 'smg', 'vwfa'};
    ROIs_coord = [-45 -24 8;-48 -44 23; -45 -51 -10];
end

nregions    = length(reg);

if nroi == 3
    lstgm=1; smg=2; vwfa=3;
elseif nroi == 4
    lstgm=1; smg=2; sts=3; vwfa=4;
end
% this code is to check how much of the variance in your VOI is explained by the timeseries :
% load('VOI_****_adapted_1.mat')
% 100*xY.s(1)/sum(xY.s)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory containing the data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nroi == 3
    AUD_dcm_matrix_3roi_design;
elseif nroi == 4
    if lexic
        AUD_dcm_matrix_4roi_design_lexic
    elseif ~lexic
        AUD_dcm_matrix_4roi_design;
    end
end

temp_a=a;
clear a

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
    
    SPM = load(fullfile(path_to_stats,dcm_folder,'SPM.mat'));
    SPM = SPM.SPM;
    
    if nroi == 3
        reg         = {'lpstg', 'smg', 'vwfa'};
    elseif nroi == 4
        reg         = {'lstgm', 'smg', 'sts', 'vwfa'};
    end
    
    nregions    = length(reg);
    
    if nroi == 3
        lstgm=1; smg=2; vwfa=3;
    elseif nroi == 4
        lstgm=1; smg=2; sts=3; vwfa=4;
    end
    
    if nroi == 3
        f = {
        fullfile(path_to_stats,sprintf('VOI_lSTGm_-45_-24_8_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_SMG_-48_-44_23_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_VWFA_-45_-51_-10_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        };
    elseif nroi == 4
        f = {
        fullfile(path_to_stats,sprintf('VOI_lSTGm_-45_-24_8_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_SMG_-48_-44_23_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_lSTS_-65_-41_6_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        fullfile(path_to_stats,sprintf('VOI_VWFA_-45_-51_-10_%dmm_sph_DCM_ROI_adapted_to_aud_con16_adj_eoi5_1.mat', sphere_radius));
        };
    end
    
    cd (path_to_stats)
    for r = 1 : nregions
       XY = load(f{r});
       xY(r) = XY.xY;
    end
    
    if lexic
        res_path = fullfile(path_to_stats, dcm_folder, 'lex_cond');
    else
        res_path = fullfile(path_to_stats, dcm_folder, 'no_lex_cond');
    end
    
    if nroi == 3
        res_path = fullfile(res_path, sprintf('all_3_rois_models_%dmm', sphere_radius));
    elseif nroi == 4
        res_path = fullfile(res_path, sprintf('all_4_rois_models_%dmm', sphere_radius));
    end
    
    if ~isdir(res_path)
        mkdir(res_path)
    end
    
    cd(res_path);
    
    clear a
    
    %% common parameters for all models;
    
    if nroi == 3
        include = [1 1 0];
    elseif nroi == 4
        if lexic
            include = [1 1 1];
        elseif ~lexic
            include = [1 1 0];
        end
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
        a=temp_a{1,i}(:,:,1);
        modul_models = length(find(~cellfun('isempty',temp_a(1:end,i)))==1);
        for j=1:modul_models
            if i>1
                name_count = length(find(~cellfun(@isempty,temp_a(1:end,1:i-1)))) + length(find(~cellfun(@isempty,temp_a(1:j,i))));
            else
                name_count = length(find(~cellfun(@isempty,temp_a(1:j,i))));
            end
            
            b(:,:,speech)=temp_a{j,i}(:,:,2);
            if nroi == 4 & lexic
                b(:,:,words)=temp_a{j,i}(:,:,3);
            end
            
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