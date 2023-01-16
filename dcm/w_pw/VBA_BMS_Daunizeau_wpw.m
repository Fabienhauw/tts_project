clear;clc;
addpath('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/scripts_finaux/dcm/w_pw');
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB/spm12')) 

i=0;
D = '/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI';
S = dir(D);
mask = ismember({S.name}, {'.', '..','FLC (copy)', 'FLC', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
counter = 0;

for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end
for k = 1:22 %parfor 1:numel(S)
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    dcm_res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models/fam_comparison');
    cd(dcm_res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L(:,k) = BMS.DCM.ffx.family.like(:);
end

% L: Kxn array of log-model evidences (K models; n subjects)
[posterior,out] = VBA_groupBMC(L) ;

clear L
for k = 1:22 %parfor 1:numel(S)
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    dcm_res_path = fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models/struc_9');
    cd(dcm_res_path)
    BMS = load('BMS.mat');
    BMS = BMS.BMS;
    L(:,k) = BMS.DCM.ffx.F(:);
end

[posterior,out] = VBA_groupBMC(L) ;

%% Bayesian model average for structural connectivity
post_prob = posterior.r';
nreg = 4;
mean_conn(1:numel(path_to_scan)) = {zeros(nreg,nreg)};
mean_modul(1:numel(path_to_scan)) = {zeros(nreg,nreg)};
for k = 1:numel(path_to_scan)
    path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
    cd(fullfile(path_to_stats, 'dcm_minus_Sujet08/all_models/struc_9'));
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

cd(fullfile('/network/lustre/dtlake01/cohen/data/Fabien_official/flc_vs_controls/fmri/FinalMRI/FLC/Aud/stats/run/dcm_model_param_modul_endroit/dcm_minus_Sujet08/all_models/struct_9'));
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


for x = 1:size(connex,1)
    for y =1:size(connex,2)
        for param = 1:numel(path_to_scan)
            conn_distrib{x,y}(param) = mean_conn{1,param}(x,y);
            modul_distrib{x,y}(param) = mean_modul{1,param}(x,y,2);
        end
    end
end


%% Plot for connections

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
            [hypoth, p, t, df] = crawford_ttest(mean_conn_flc(x,y), conn_distrib{x,y});
            p_values(x,y)=p;
            title(sprintf('Connexion from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);


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
            [hypoth, p, t, df] = crawford_ttest(mean_modul_flc(x,y,2), modul_distrib{x,y});
            p_values(x,y)=p;
            title(sprintf('Modulation from %s to %s', reg{y}, reg{x}));
        end
    end
end
p_values_res = reshape(p_values, 1, numel(p_values)); mask = find(p_values_res == 0); p_values_res(mask) = [];
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values_res);