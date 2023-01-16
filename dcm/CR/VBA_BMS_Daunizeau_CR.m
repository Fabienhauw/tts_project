% Comparison of family model evidences

clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/CR');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC', 'FLC (copy)', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
path_to_all_stats = {};
counter = 0;

for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/rosso_controls';
S = dir(D);
mask = ismember({S.name}, {'.', '..', 'README'});
S(mask) = [];

nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs               = [subjs;S(k).name];
    path_to_subj        = fullfile(D, S(k).name);
    path_to_scan        = [path_to_scan;path_to_subj];
    tmp_path            = fullfile(path_to_subj, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    path_to_all_stats   = [path_to_all_stats; tmp_path];
end

for k = 1:35
    path_to_stats = path_to_all_stats{k};
    dcm_res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models/less_modul/fam_comparison');
    cd(dcm_res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L(:,k) = BMS.DCM.ffx.family.like(:);
end

% L: Kxn array of log-model evidences (K models; n subjects)
[posterior,out] = VBA_groupBMC(L) ;

clear L
for k = 1:35 %parfor 1:numel(S)
    path_to_stats=path_to_all_stats{k};
%     try rmdir(fullfile(path_to_stats,'dcm_minus_Sujet08/all_models/less_modul/struc_9'),'s'), catch, end;
    dcm_res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models/less_modul/struc_9');
    cd(dcm_res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L(:,k) = BMS.DCM.ffx.F(:);
end

[posterior,out] = VBA_groupBMC(L) ;

%% Bayesian model average for structural connectivity
post_prob = posterior.r';
nreg = 4;
mean_conn(1:numel(path_to_all_stats)) = {zeros(nreg,nreg)};
mean_modul(1:numel(path_to_all_stats)) = {zeros(nreg,nreg)};
for k = 1:numel(path_to_all_stats)
    cd(fullfile(path_to_all_stats{k}, 'dcm_minus_Sujet08/all_models/less_modul/struc_9'));
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
        mean_conn{k}    = [mean_conn{k}+wconnex]; % this is the sum for each subj for each connection parameter in each model, weighted by the posterior probability of the model
        mean_modul{k}   = [mean_modul{k}+wmodul]; % same for modulation
    end
end

cd(fullfile('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI/FLC/Aud/stats/chaperon/dcm_model_param_modul_endroit/dcm_minus_Sujet08/all_models/less_modul/struct_9'));
tmp_post_prob = load('BMS.mat');
tmp_post_prob = tmp_post_prob.BMS.DCM.ffx.model.post;
subj = load('model_space.mat'); subj = subj.subj;
models(:) = {subj.sess.model(:).fname};
mean_conn_flc = zeros(nreg,nreg);
mean_modul_flc = zeros(nreg,nreg);
for dcm_i = 1:length(models)
    tmp_DCM_flc  = load(models{dcm_i});
    tmp_DCM_flc     = tmp_DCM_flc.DCM;
    connex_flc      = tmp_DCM_flc.Ep.A;
    modul_flc       = tmp_DCM_flc.Ep.B;
    wconnex_flc     = tmp_post_prob(dcm_i)*connex_flc;
    wmodul_flc      = tmp_post_prob(dcm_i)*modul_flc;
    mean_conn_flc    = [mean_conn_flc+wconnex_flc];
    mean_modul_flc   = [mean_modul_flc+wmodul_flc];
end

vector_cov1 = [28;26;30;26;29;29;29;28;36;29;24;20;29;30;28;28;20;20;23;27;27;...
    61;69;61;64;72;67;75;60;68;61;69;64;52;51;69];

for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(path_to_all_stats)
            conn_distrib{x,y}(param) = mean_conn{1,param}(x,y);
            modul_distrib{x,y}(param) = mean_modul{1,param}(x,y,2);
        end
    end
end

[cb] = cbrewer('qual', 'Set1', 16, 'pchip'); reg = {'STGm', 'SMG', 'MTG', 'VWFA'};
count = 0;
for x = 1:size(connex,1)
    for y =1:size(connex,2)
        if any(conn_distrib{x,y})~=0
            count = count + 1;
            f(count) = figure;
            fig{count} = raincloud_plot(conn_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
                'box_col_match', 1);
            line = xline(mean_conn_flc(x,y), '-.r');
            set(line,'linewidth',3)
            
            % stats
            fitlm_conn = fitlm(vector_cov1,[conn_distrib{x,y} mean_conn_flc(x,y)]);
            fitlm_conn = fitlm_conn.Residuals.Raw;
            conn_distrib_fit{x,y} = fitlm_conn(1:end-1)'; mean_conn_flc_fit(x,y) = fitlm_conn(end);

            [hypoth, p, t, df] = crawford_ttest(mean_conn_flc_fit(x,y), conn_distrib_fit{x,y});
            p_values(x,y)=p;
            title(sprintf('Connexion from %s to %s', reg{y}, reg{x}));
            
            sd(x,y) = std(conn_distrib{x,y});
            nb_sd(x,y) = (mean_conn_flc(x,y)-mean(conn_distrib{x,y}))/sd(x,y);
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

return

%% Plot for modulations
clear p_values
count = 0;
for x = 1:size(modul(:,:,2),2)
    for y =1:size(modul(:,:,2),2)
        if any(modul_distrib{x,y})~=0
            count = count + 1;
            f(count) = figure;
            fig{count} = raincloud_plot(modul_distrib{x,y}, 'box_on', 1, 'color', cb(count,:), 'alpha', 0.5,...
                'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
                'box_col_match', 1);
            line = xline(mean_modul_flc(x,y,2), '-.r');
            set(line,'linewidth',3)
            title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
            
            % stats
            fitlm_mod = fitlm(vector_cov1,[modul_distrib{x,y} mean_modul_flc(x,y,2)]);
            fitlm_mod = fitlm_mod.Residuals.Raw;
            modul_distrib_fit{x,y} = fitlm_mod(1:end-1); mean_modul_flc_fit(x,y) = fitlm_mod(end);
            
            [hypoth, p, t, df] = crawford_ttest(mean_modul_flc_fit(x,y), modul_distrib_fit{x,y});
            p_values(x,y)=p;
            
            sd(x,y) = std(modul_distrib{x,y});
            nb_sd(x,y)=(mean_modul_flc(x,y,2)-mean(modul_distrib{x,y}))/sd(x,y);
        end
    end
end

p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);

