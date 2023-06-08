clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')) 

nroi = 4;
lexic = 1;
fam_comp = 1;
redo_comp = 0;
ang_gyr = 1;
model_kind = 3;
if model_kind == 1
    dcm_folder = 'dcm_model_param_modul_speech_baseline';
elseif model_kind == 2
    dcm_folder = 'dcm_model_param_modul_sent_scramble';
elseif model_kind == 3
    dcm_folder = 'dcm_model_param_modul';
end

sphere_radius = 4;
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
if nroi == 3
    AUD_dcm_matrix_3roi_design;
elseif nroi == 4 & ~lexic
    AUD_dcm_matrix_4roi_design;
elseif nroi == 4 & lexic
    AUD_dcm_matrix_4roi_design_lexic;
end
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
    if lexic & ~ang_gyr
        tmp_final_path = fullfile(tmp_path, dcm_folder, 'lex_cond');
    elseif lexic & ang_gyr
        tmp_final_path = fullfile(tmp_path, dcm_folder, 'ang_gyr');
    elseif ~lexic & ~ang_gyr
        tmp_final_path = fullfile(tmp_path, dcm_folder, 'no_lex_cond');
    end
    if nroi == 3
        tmp_final_path = fullfile(tmp_final_path,sprintf('all_3_rois_models_%dmm', sphere_radius));
    elseif nroi == 4
        tmp_final_path = fullfile(tmp_final_path,sprintf('all_4_rois_models_%dmm', sphere_radius));
    end
    final_path          = [final_path ;tmp_final_path];
end


for k = 1 : numel(S)
    fprintf('%s \n', S(k).name)
    res_path = final_path{k};

    cd(res_path);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         MODEL COMPARISON           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if redo_comp
        try
            delete('BMS.mat');
        catch
        end
    end
    all_DCM = dir('DCM*');
    path_to_dcm = {};
    for dcm=1:length(all_DCM)
        new_dcm      = fullfile(res_path,all_DCM(dcm).name);
        path_to_dcm  = [path_to_dcm;new_dcm];
    end
    
    %to compare all DCM models;
    clear matlabbatch
    fprintf('All models comparison \n')

    matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path};
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm;
    matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
    matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
    matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
    matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
    matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
    
    if ~exist(fullfile(res_path, 'BMS.mat')) | redo_comp
        spm_jobman('run', matlabbatch)
    end
    clear matlabbatch
    
    
    fprintf('All families comparison \n')
    if size(temp_a,2)>1 & fam_comp
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
        
        cd(fullfile(res_path,'fam_comparison'))
        if redo_comp
            try
                delete('BMS.mat');
            catch
            end
        end
        
        matlabbatch{1}.spm.dcm.bms.inference.dir = {fullfile(res_path,'fam_comparison')};
        matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm;
        matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
        matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
        matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
        matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {fullfile(res_path,'family.mat')};
        matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
        matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
        
        if ~exist(fullfile(res_path, 'fam_comparison/BMS.mat')) | redo_comp
            spm_jobman('run', matlabbatch)
        end
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
        cd(res_path)
        DCM_winning_fam = dir(sprintf('*%s*.mat',win_fam));
        path_to_dcm_win_fam = {};
        for dcm=1:length(DCM_winning_fam)
            new_dcm      = fullfile(res_path,DCM_winning_fam(dcm).name);
            path_to_dcm_win_fam= [path_to_dcm_win_fam;new_dcm];
        end
        
        %new DCM comparison:
        fprintf('All models among best family comparison \n')
        cd(res_path_fam)
        if redo_comp
            try
                delete('BMS.mat');
            end
        end
        matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam};
        matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_win_fam;
        matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
        matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
        matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
        matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
        matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
        matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
        
        if ~exist(fullfile(res_path_fam, 'BMS.mat')) | redo_comp
            spm_jobman('run', matlabbatch)
        end
        clear matlabbatch

    end
    
    %% to compare systematically all models within the best group model (ie
    % struct 9)
    
    path_to_dcm_fam_9 = {};
    cd(res_path)
    %identify the DCM within the 9th family:
    DCM_fam_9 = dir(sprintf('*struct_9*.mat'));
    for dcm=1:length(DCM_fam_9)
        new_dcm             = fullfile(res_path,DCM_fam_9(dcm).name);
        path_to_dcm_fam_9   = [path_to_dcm_fam_9;new_dcm];
    end
    
    res_path_fam_9 = fullfile(res_path, 'struct_9');
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
    if ~exist(fullfile(res_path_fam_9, 'BMS.mat')) | redo_comp
        spm_jobman('run', matlabbatch)
    end
    clear matlabbatch
end

%% to make a .mat with all info for each subject about the 3 previous comparisons

all_results(numel(S)).name                  = [];
all_results(numel(S)).bestmodel             = [];
all_results(numel(S)).bestfam               = [];
all_results(numel(S)).bestmodel_inbestfam   = [];

for k = 1 : numel(S)
    fprintf('%s \n', S(k).name)
    all_results(k).name = S(k).name;
    res_path = final_path{k};
    
    
    cd(res_path);
    BMS = load('BMS.mat'); models = load('model_space.mat');
    BMS = BMS.BMS;
    [val, idx] = max(BMS.DCM.ffx.F);
    path_bestmodel = strsplit(models.subj.sess.model(idx).fname, '/'); name_bestmodel = path_bestmodel{end};
    all_results(k).bestmodel = name_bestmodel;
    
    BMS = load(fullfile(res_path,'fam_comparison/BMS.mat'));
    BMS = BMS.BMS;
    [val, idx] = max(BMS.DCM.ffx.family.post);
    win_fam         = sprintf('struct_%d',idx);
    all_results(k).bestfam = win_fam;
    
    BMS = load(fullfile(res_path, win_fam, 'BMS.mat')); models = load(fullfile(res_path, win_fam, 'model_space.mat'));
    BMS = BMS.BMS;
    [val, idx] = max(BMS.DCM.ffx.F);
    path_bestmodel = strsplit(models.subj.sess.model(idx).fname, '/'); name_bestmodel = path_bestmodel{end};
    all_results(k).bestmodel_inbestfam = name_bestmodel;

end