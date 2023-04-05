clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')) 

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

counter = 0;
clear a
AUD_dcm_matrix_3roi_design;
% AUD_dcm_matrix_4roi_design;
temp_a=a;
clear a

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


for k = 1:numel(S)
    res_path=fullfile(path_to_all_stats{k}, 'dcm_model_param_modul/all_3_rois_models_4mm');
%     res_path=fullfile(path_to_all_stats{k}, 'dcm_model_param_modul/all_4_rois_models_4mm');

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
    
    if size(temp_a,2)>1
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
    end
    
%     %% to compare systematically all models within the best group model (ie
%     % struct 9)
%     
%     path_to_dcm_fam_9 = {};
%     %identify the DCM within the 9th family:
%     DCM_fam_9 = dir(sprintf('*struct_9*.mat'));
%     for dcm=1:length(DCM_fam_9)
%         new_dcm             = fullfile(res_path,DCM_fam_9(dcm).name);
%         path_to_dcm_fam_9   = [path_to_dcm_fam_9;new_dcm];
%     end
%     
%     res_path_fam_9 = fullfile(res_path, 'struc_9');
%     if ~isdir(res_path_fam_9)
%         mkdir(res_path_fam_9)
%     end
% 
%     %new DCM comparison:
%     matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam_9};
%     matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_fam_9;
%     matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
%     matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
%     matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
%     matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
%     matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
%     matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
%     
%     spm_jobman('run', matlabbatch)
%     clear matlabbatch
end