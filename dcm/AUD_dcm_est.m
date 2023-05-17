clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nroi = 4;
lexic = 1;
redo =0;
model_kind = 2;
if model_kind == 1
    dcm_folder = 'dcm_model_param_modul_speech_baseline';
elseif model_kind == 2
    dcm_folder = 'dcm_model_param_modul_sent_scramble';
elseif model_kind == 3
    dcm_folder = 'dcm_model_param_modul';
end

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
counter = 0;

for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/loc/stats_s5_without_resting');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

for k = 1 : numel(S)
    path_to_stats = path_to_all_stats{k};
    
    if lexic
        res_path = fullfile(path_to_stats, dcm_folder, 'lex_cond');
    else
        res_path = fullfile(path_to_stats, dcm_folder, 'no_lex_cond');
    end
    
    if nroi == 3
        res_path = fullfile(res_path, 'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path, 'all_4_rois_models_4mm');
    end
    
    if ~isdir(res_path)
        mkdir(res_path)
    end
    cd(res_path)

%     all_DCM_1 = dir('DCM_struct_2*'); all_DCM_2 = dir('DCM_struct_3*');
%     all_DCM = [all_DCM_1;all_DCM_2];
        all_DCM = dir('DCM_*');
    path_to_dcm = {};
    all_DCM_non_estimated = {};
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    for dcm=1:length(all_DCM)
        try
            tmp_dcm = load(all_DCM(dcm).name);
        catch
            disp(['still impossible to load: ' fullfile(res_path,all_DCM(dcm).name) ', try to estimate it manually...'])
        end
        
        if ~isfield(tmp_dcm, 'F')
            all_DCM_non_estimated = [all_DCM_non_estimated;fullfile(res_path,all_DCM(dcm).name)];
%             disp(['non estimated: ' fullfile(res_path,all_DCM(dcm).name)])
        end
    end
    
    clear a
    clear b
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL ESTIMATION           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %DCM Estimation
    if redo
        for dcm_mod = 1:length(path_to_dcm)
            counter = counter + 1;
            matlabbatch{counter}.spm.dcm.estimate.dcms.subj.dcmmat = path_to_dcm(dcm_mod);
            matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
        end
    elseif ~redo
        for dcm_mod = 1:length(all_DCM_non_estimated)
            counter = counter + 1;
            matlabbatch{counter}.spm.dcm.estimate.dcms.subj.dcmmat = all_DCM_non_estimated(dcm_mod);
            matlabbatch{counter}.spm.dcm.estimate.output.separate = struct([]);
        end
    end
    %%%%%%%% this line below to comment to avoid re estimating
    %%%%%%%% models
    
end

%%
cd('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm')

par.run = 0;
par.sge = 1;
par.sge_queu = 'normal,bigmem';
par.pct = 1;
par.walltime = '00:30:00';
par.jobname  = 'dcm_estimate';
%%%%%%%% this line below to comment to avoid re estimating
%%%%%%%% models

job_ending_rountines(matlabbatch, [], par);
