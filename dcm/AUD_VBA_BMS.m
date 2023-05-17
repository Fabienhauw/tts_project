clear;clc;
addpath('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm');
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nroi = 4;
lexic = 1;

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

%% here we estimate posterior probability at the group level for each model based on probability for each individuals models comparison
clear L1 L2

for k = 1:numel(S)/2
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path)
    BMS = load('BMS.mat'); BMS = BMS.BMS;
    L1(:,k) = BMS.DCM.ffx.F(:);
end

[posterior1,out1] = VBA_groupBMC(L1) ;

for k = numel(S)/2+1:numel(S) %parfor 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L2(:,k-(numel(S)/2)) = BMS.DCM.ffx.F(:);
end
[posterior2,out2] = VBA_groupBMC(L2) ;


%--------------------------------------------------------------------------


% Bayesian model average for structural connectivity
% for group 1
post_prob = posterior1.r'; %% quoi utiliser ? posterior.r ou .a???? Revoir avec Jean
mean_conn1(1:numel(S)/2) = {zeros(nroi,nroi)};
mean_modul1(1:numel(S)/2) = {zeros(nroi,nroi)};
clear conn_distrib1 modul1_distrib1 modul2_distrib1 conn_distrib2 modul1_distrib2 modul2_distrib2
for k = 1:numel(S)/2
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path);
    tmp_post_prob = post_prob(k,:);
    subj = load('model_space.mat'); subj = subj.subj;
    models(:) = {subj.sess.model(:).fname};
    for dcm_i = 1:length(models)
        tmp_DCM     = load(models{dcm_i});
        tmp_DCM     = tmp_DCM.DCM;
        connex      = tmp_DCM.Ep.A;
        modul       = tmp_DCM.Ep.B;
        wconnex     = tmp_post_prob(dcm_i)*connex;
        wmodul      = tmp_post_prob(dcm_i)*modul;
        mean_conn1{k}    = [mean_conn1{k}+wconnex];
        mean_modul1{k}   = [mean_modul1{k}+wmodul];
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)/2
            conn_distrib1{x,y}(param) = mean_conn1{1,param}(x,y);
            modul1_distrib1{x,y}(param) = mean_modul1{1,param}(x,y,2);
            modul2_distrib1{x,y}(param) = mean_modul1{1,param}(x,y,3);
        end
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)/2
            mconn_distrib1{x,y} = mean(conn_distrib1{x,y});
            mmodul1_distrib1{x,y} = mean(modul1_distrib1{x,y});
            mmodul2_distrib1{x,y} = mean(modul2_distrib1{x,y});
        end
    end
end

% for group 2
post_prob = posterior2.r';
mean_conn2(1:numel(S)/2) = {zeros(nroi,nroi)};
mean_modul2(1:numel(S)/2) = {zeros(nroi,nroi)};
for k = numel(S)/2+1:numel(S)
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path);
    tmp_post_prob = post_prob(k-numel(S)/2,:);
    subj = load('model_space.mat'); subj = subj.subj;
    models(:) = {subj.sess.model(:).fname};
    for dcm_i = 1:length(models)
        tmp_DCM     = load(models{dcm_i});
        tmp_DCM     = tmp_DCM.DCM;
        connex      = tmp_DCM.Ep.A;
        modul       = tmp_DCM.Ep.B;
        wconnex     = tmp_post_prob(dcm_i)*connex;
        wmodul      = tmp_post_prob(dcm_i)*modul;
        mean_conn2{k-numel(S)/2}    = [mean_conn2{k-numel(S)/2}+wconnex];
        mean_modul2{k-numel(S)/2}   = [mean_modul2{k-numel(S)/2}+wmodul];
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)/2
            conn_distrib2{x,y}(param) = mean_conn2{1,param}(x,y);
            modul1_distrib2{x,y}(param) = mean_modul2{1,param}(x,y,2);
            modul2_distrib2{x,y}(param) = mean_modul2{1,param}(x,y,3);
        end
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)/2
            mconn_distrib2{x,y} = mean(conn_distrib2{x,y});
            mmodul1_distrib2{x,y} = mean(modul1_distrib2{x,y});
            mmodul2_distrib2{x,y} = mean(modul2_distrib2{x,y});
        end
    end
end


% Plot for connections

% [cb] = cbrewer('qual', 'Set1', 16, 'pchip'); reg = {'STGm', 'SMG', 'MTG', 'VWFA'};
count = 0;
for x = 1:size(connex,1)
    for y =1:size(connex,2)
%         if any(conn_distrib1{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(conn_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_conn_flc(x,y), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(conn_distrib1{x,y}, conn_distrib2{x,y});
            p_values(x,y)=p;
%             title(sprintf('Connexion from %s to %s', reg{y}, reg{x}));
%         end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);


% Plot for modulations 1

clear p_values
count = 0;
for x = 1:size(modul(:,:,2),2)
    for y =1:size(modul(:,:,2),2)
        if any(modul1_distrib1{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(modul_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_modul_flc(x,y,2), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(modul1_distrib1{x,y}, modul1_distrib2{x,y});
            p_values(x,y)=p;
%             title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

% Plot for modulations 2

clear p_values
count = 0;
for x = 1:size(modul(:,:,2),2)
    for y =1:size(modul(:,:,2),2)
        if any(modul2_distrib1{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(modul_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_modul_flc(x,y,2), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(modul2_distrib1{x,y}, modul2_distrib2{x,y});
            p_values(x,y)=p;
%             title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% New estimation if estimation of all models probability is made with all subjects as one group
clear L

for k = 1:numel(S) %parfor 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L(:,k) = BMS.DCM.ffx.F(:);
end

[posterior,out] = VBA_groupBMC(L) ;

%--------------------------------------------------------------------------
% Bayesian model average for structural connectivity
% for group 1
post_prob = posterior.r';
mean_conn(1:numel(S)) = {zeros(nroi,nroi)};
mean_modul(1:numel(S)) = {zeros(nroi,nroi)};
clear conn_distrib modul_distrib
for k = 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    
    cd(res_path);
    tmp_post_prob = post_prob(k,:);
    subj = load('model_space.mat'); subj = subj.subj;
    models(:) = {subj.sess.model(:).fname};
    for dcm_i = 1:length(models)
        tmp_DCM     = load(models{dcm_i});
        tmp_DCM     = tmp_DCM.DCM;
        connex      = tmp_DCM.Ep.A;
        modul       = tmp_DCM.Ep.B;
        wconnex     = tmp_post_prob(dcm_i)*connex;
        wmodul      = tmp_post_prob(dcm_i)*modul;
        mean_conn{k}    = [mean_conn{k}+wconnex];
        mean_modul{k}   = [mean_modul{k}+wmodul];
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)
            conn_distrib{x,y}(param) = mean_conn{1,param}(x,y);
            modul1_distrib{x,y}(param) = mean_modul{1,param}(x,y,2);
            modul2_distrib{x,y}(param) = mean_modul{1,param}(x,y,3);
        end
    end
end

%--------------------------------------------------------------------------
% Plot for connections

% [cb] = cbrewer('qual', 'Set1', 16, 'pchip'); reg = {'STGm', 'SMG', 'MTG', 'VWFA'};
count = 0;
for x = 1:size(connex,1)
    for y =1:size(connex,2)
        if any(conn_distrib{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(conn_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_conn_flc(x,y), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(conn_distrib{x,y}(1:17), conn_distrib{x,y}(18:end));
            p_values(x,y)=p;
%             title(sprintf('Connexion from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);


% Plot for modulations 1

clear p_values
count = 0;
for x = 1:size(modul(:,:,2),2)
    for y =1:size(modul(:,:,2),2)
        if any(modul1_distrib{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(modul_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_modul_flc(x,y,2), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(modul1_distrib{x,y}(1:17), modul1_distrib{x,y}(18:end));
            p_values(x,y)=p;
%             title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

% Plot for modulations 2

clear p_values
count = 0;
for x = 1:size(modul(:,:,2),2)
    for y = 1 : size(modul(:,:,2),2)
        if any(modul2_distrib{x,y})~=0
            count = count + 1;
%             f(count) = figure;
%             fig{count} = raincloud_plot(modul_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
%                 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%                 'box_col_match', 1);
%             line = xline(mean_modul_flc(x,y,2), '-.r');
%             set(line,'linewidth',3)
            [hypoth, p, t, df] = ttest2(modul2_distrib{x,y}(1:17), modul2_distrib{x,y}(18:end));
            p_values(x,y)=p;
%             title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);



%% make the sum of each weighted parameter by the model probability in each subject
% here we take the BMS and model space that was estimated across all models
% for each subject, and extracting weighted parameters based on that.

clear L
mean_conn(1:numel(S)) = {zeros(nroi,nroi)};
mean_modul(1:numel(S)) = {zeros(nroi,nroi)};
clear conn_distrib modul_distrib

for k = 1:numel(S) %parfor 1:numel(S)
    path_to_stats = path_to_all_stats{k};
    if lexic
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/lex_cond');
    else
        res_path = fullfile(path_to_stats, 'dcm_model_param_modul/no_lex_cond');
    end
    if nroi == 3
        res_path = fullfile(res_path,'all_3_rois_models_4mm');
    elseif nroi == 4
        res_path = fullfile(res_path,'all_4_rois_models_4mm');
    end
    cd(res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    posterior(:,k) = BMS.DCM.ffx.model.post(:)';
    
    tmp_post_prob = posterior(:,k);
    subj = load('model_space.mat'); subj = subj.subj;
    models(:) = {subj.sess.model(:).fname};
    for dcm_i = 1:length(models)
        tmp_DCM     = load(models{dcm_i});
        tmp_DCM     = tmp_DCM.DCM;
        connex      = tmp_DCM.Ep.A;
        modul       = tmp_DCM.Ep.B;
        wconnex     = tmp_post_prob(dcm_i)*connex;
        wmodul      = tmp_post_prob(dcm_i)*modul;
        mean_conn{k}    = [mean_conn{k}+wconnex];
        mean_modul{k}   = [mean_modul{k}+wmodul];
    end
end

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(S)
            conn_distrib{x,y}(param) = mean_conn{1,param}(x,y);
            modul1_distrib{x,y}(param) = mean_modul{1,param}(x,y,2);
            modul2_distrib{x,y}(param) = mean_modul{1,param}(x,y,3);
        end
    end
end

%% stats for connections and modulations

count = 0;
for x = 1:size(connex,1)
    for y =1:size(connex,2)
        if any(conn_distrib{x,y})~=0
            count = count + 1;
            [hypoth, p, t, df] = ttest2(conn_distrib{x,y}(1:17), conn_distrib{x,y}(18:end));
            p_values1(x,y)=p;
            [hypoth, p, t, df] = ttest2(modul1_distrib{x,y}(1:17), modul1_distrib{x,y}(18:end));
            p_values2(x,y)=p;
            [hypoth, p, t, df] = ttest2(modul2_distrib{x,y}(1:17), modul2_distrib{x,y}(18:end));
            p_values3(x,y)=p;
        end
    end
end
p_values_res = reshape(p_values1, 1, numel(p_values1)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

p_values_res = reshape(p_values2, 1, numel(p_values2)); mask = find(p_values_res == 0 | isnan(p_values_res)); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

p_values_res = reshape(p_values3, 1, numel(p_values3)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);