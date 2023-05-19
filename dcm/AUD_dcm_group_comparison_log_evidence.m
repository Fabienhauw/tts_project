clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/VBA-toolbox'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nroi = 4;
lexic = 1;
% redo_comp = 1;
model_kind = 1;
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

%% First part is to made the matrix L1 for group 1 and L2 for group 2, that will includ the evidence for each family for each participant
clear L1 L2 a

if nroi == 3
    AUD_dcm_matrix_3roi_design;
elseif nroi == 4
    AUD_dcm_matrix_4roi_design;
end

temp_a=a;

final_fam_evid = [];
for k = 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    disp(S(k).name)
    if lexic
        res_path = fullfile(path_to_stats, dcm_folder, 'lex_cond');
    else
        res_path = fullfile(path_to_stats, dcm_folder, 'no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path, '/all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path, 'all_4_rois_models_4mm');
    end
    cd(res_path)
    BMS = load('BMS.mat'); BMS = BMS.BMS;
    
    models = load('model_space.mat'); 
    all_models = {};
    for tmp_model = 1 : size(models.subj.sess.model,2) 
        all_models = [all_models;models.subj.sess.model(tmp_model).fname];
    end
    
    for tmp_fam = 1 : size(temp_a,2)
        tmp_mask = ~cellfun(@isempty,(regexp(all_models,sprintf('struct_%d', tmp_fam))));
        fam_models = all_models(tmp_mask);
        fam_evid = BMS.DCM.ffx.F(tmp_mask);
        norm_factor = max(fam_evid);
        
        fam_prob = 1/size(fam_models,1);
        sum_model_proba = 0;
        for tmp_model = 1 : size(fam_models,1)
            model_sub_proba = fam_prob*exp(fam_evid(tmp_model) - norm_factor); %proba of each model * exp(evidence of this model)
            sum_model_proba = sum_model_proba + model_sub_proba;
        end
        
        final_fam_evid(tmp_fam,k) = norm_factor + log(sum_model_proba);
    end
end

final_fam_evid_g1 = final_fam_evid(:, 1:17);
final_fam_evid_g2 = final_fam_evid(:, 18:end);

[posterior,out] = VBA_groupBMC(final_fam_evid) ;
Fe = out.F(end);

[posterior1, out1] = VBA_groupBMC(final_fam_evid_g1) ;
[posterior2, out2] = VBA_groupBMC(final_fam_evid_g2) ;
Fd = out1.F(end) + out2.F(end) ;

p = 1/(1+exp(Fd-Fe)); % p> 0.5, it is most likely that both groups have the same best family;


%% if same families, then do the same analysis on models scale;
if p > 0.5
    clear L
    [val, idx] = max(out.Ef);
    for k = 1:numel(S)
        path_to_stats = path_to_all_stats{k};
        disp(S(k).name)
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
        cd(res_path)
%         idx = 9;
        win_fam         = sprintf('struct_%d',idx);
        res_path_fam = fullfile(res_path, win_fam);
        if ~isdir(res_path_fam)
            mkdir(res_path_fam)
        end
        
        cd(res_path_fam)
%         if redo_comp
%             delete('BMS.mat');
%         end
        
        cd(res_path)

        %identify the DCM within the winning family:
        DCM_winning_fam = dir(sprintf('*%s*.mat',win_fam));
        path_to_dcm_win_fam = {};
        for dcm=1:length(DCM_winning_fam)
            new_dcm      = fullfile(res_path,DCM_winning_fam(dcm).name);
            path_to_dcm_win_fam= [path_to_dcm_win_fam;new_dcm];
        end
        
        % estimate the evidence for each model for each participant in the
        % best family
        fprintf('All models among best family comparison \n')
        matlabbatch{1}.spm.dcm.bms.inference.dir = {res_path_fam};
        matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = path_to_dcm_win_fam;
        matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
        matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
        matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
        matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
        matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
        matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
        
%         spm_jobman('run', matlabbatch)
        clear matlabbatch
        
        cd(res_path_fam)
        BMS = load('BMS.mat'); BMS = BMS.BMS;
        L(:,k) = BMS.DCM.ffx.F(:);
        
    end
    
    L1 = L(:, 1:17);
    L2 = L(:, 18:end);
    
    [posterior,out] = VBA_groupBMC(L);
    Fe = out.F(end);
    
    [posterior1, out1] = VBA_groupBMC(L1) ;
    [posterior2, out2] = VBA_groupBMC(L2) ;
    Fd = out1.F(end) + out2.F(end) ;
    
    p = 1/(1+exp(Fd-Fe)); % p> 0.5, it is most likely that both groups have the same best model;
    
end

if p >0.5
    
    [val, idx] = max(out.Ef); % here it is the fam 9, num 55;
    winning_model_split = strsplit(path_to_dcm_win_fam{idx}, '/');
    winning_model = winning_model_split{end};
end

%% Last part would be to test the mean parameters for this best model
if p > 0.5
    for k = 1:numel(S)
        path_to_stats = path_to_all_stats{k};
        disp(S(k).name)
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
        cd(res_path)
        best_model = winning_model;
        tmp_DCM     = load(best_model); tmp_DCM     = tmp_DCM.DCM;
        connex{k}   = tmp_DCM.Ep.A;
        modul{k}    = tmp_DCM.Ep.B;
    end
    
    count = 0;
    for x = 1:size(connex{1,1},1)
        for y =1:size(connex{1,1},2)
            for param = 1:numel(S)
                conn_distrib{x,y}(param) = connex{1,param}(x,y);
                modul1_distrib{x,y}(param) = modul{1,param}(x,y,2);
                modul2_distrib{x,y}(param) = modul{1,param}(x,y,3);
            end
            
            if any(conn_distrib{x,y})~=0
                count = count + 1;
                [hypoth, p, t, df] = ttest2(conn_distrib{x,y}(1:17), conn_distrib{x,y}(18:end));
                mconn_distrib1{x,y} = mean(conn_distrib{x,y}(1:17));
                mconn_distrib2{x,y} = mean(conn_distrib{x,y}(18:end));
                p_values_conn(x,y)=p;
                
                [hypoth, p, t, df] = ttest2(modul1_distrib{x,y}(1:17), modul1_distrib{x,y}(18:end));
                mmodul1_distrib1{x,y} = mean(modul1_distrib{x,y}(1:17));
                mmodul1_distrib2{x,y} = mean(modul1_distrib{x,y}(18:end));
                p_values_mod1(x,y)=p;
                
                [hypoth, p, t, df] = ttest2(modul2_distrib{x,y}(1:17), modul2_distrib{x,y}(18:end));
                mmodul2_distrib1{x,y} = mean(modul2_distrib{x,y}(1:17));
                mmodul2_distrib2{x,y} = mean(modul2_distrib{x,y}(18:end));
                p_values_mod2(x,y)=p;
            end
        end
    end
    
end

%% another option is to initially compare all models, not only the families.
clear L L1 L2
for k = 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    disp(S(k).name)
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
    cd(res_path)
    BMS = load('BMS.mat'); BMS = BMS.BMS;
    
    L(:,k) = BMS.DCM.ffx.F(:);
end

L1 = L(:, 1:17);
L2 = L(:, 18:end);

[posterior,out] = VBA_groupBMC(L);
Fe = out.F(end);

[posterior1, out1] = VBA_groupBMC(L1) ;
[posterior2, out2] = VBA_groupBMC(L2) ;
Fd = out1.F(end) + out2.F(end) ;

p = 1/(1+exp(Fd-Fe)); % p> 0.5, it is most likely that both groups have the same best family;

if p >0.5
    [val, idx] = max(out.Ef); % here it is the fam 9, num 55;
    winning_model_split = strsplit(all_models{idx}, '/');
    winning_model_across_all = winning_model_split{end};
end
