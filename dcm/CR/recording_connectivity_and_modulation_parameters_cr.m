clear;clc;
addpath('/home/fabien.hauw/MATLAB Add-Ons/spm12');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers/DCM');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers/plot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   uncomment to apply the scripts to all subjects/controls   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
repet_time = {};
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
    repet_time = [repet_time;1.3];
end

Batches.rosso.path_to_subject;
nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    repet_time = [repet_time;3];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end

%% get coordinates for all VOI used in DCM for all subjects
Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..', 'Sujet08', 'FLC_func'});
S(mask) = [];
subjs={};
nb_subj=numel(S);
path_to_scan = {};
repet_time = {};
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
    repet_time = [repet_time;1.3];
end

Batches.rosso.path_to_subject;
nb_subj=nb_subj + numel(S);
for k=1:numel(S)
    subjs=[subjs;S(k).name];
    repet_time = [repet_time;3];
    path_to_subj = fullfile(D, S(k).name);
    path_to_scan = [path_to_scan;path_to_subj];
end

for k = 1:nb_subj
    if exist(fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit');
    else
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/dcm_model_param_modul_endroit');
    end
%     path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/dcm_model');
    cd(path_to_stats)
    voi=dir('VOI_dcm*.mat');
    for v=1:length(voi)
        load(voi(v).name);
%         voi_dcm{k+1,v+1}=xY.xyz;
        if k==1
            voi_dcm{1,v+1}=xY.name;
            voi_dcm{k+1,v+1}=xY.xyz;
        else
            for dcm=1:size(voi_dcm,2)
                if strcmp(xY.name,voi_dcm{1,dcm})
                    voi_dcm{k+1,dcm}=xY.xyz;
                end
            end
        end
        if v==1
            voi_dcm{k+1,1}=subjs(k);
        end
    end
end

for k=1:nb_subj
    if sum(cellfun(@isempty,voi_dcm(k+1,2:end)))>0
        mask(k)=1;
    else
        mask(k)=0;
    end
end

path_to_scan(mask) = '';
new_nsubj = length (path_to_scan);

for k=1:new_nsubj
    if exist(fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc/dcm_model_param_modul_endroit/dcm/all_models');
    elseif isempty(regexp(path_to_scan{k}, 'FLC'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/dcm_model_param_modul_endroit/dcm/all_models');
    else
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/chaperon/dcm_model_param_modul_endroit/dcm_minus_Sujet08/all_models');
    end
    cd(path_to_stats)
    dcm = dir('DCM*struct_9_num_49*');
    load(dcm.name)

    connect{k} = DCM.Ep.A;
    modul{k} = DCM.Ep.B;
end

vector_cov1 = [69;28;26;30;26;29;29;29;28;36;29;24;20;29;30;28;28;20;20;23;27;27;...
    61;69;61;64;72;67;75;60;68;61;69;64;52;51];


modulation_1 = [];
modulation_2 = [];
conn_1 = [];
conn_2 = [];
conn_3 = [];
conn_4 = [];
conn_5 = [];

for j=1:new_nsubj
    modulation_1(j) = modul{1,j}(2,1,2);
    modulation_2(j) = modul{1,j}(3,1,2);
    conn_1(j)       = connect{1,j}(2,1);
    conn_2(j)       = connect{1,j}(3,1);
    conn_3(j)       = connect{1,j}(3,2);
%     conn_4(j)       = connect{1,j}(4,2); % this connexion is not present
%     in the new best model;
    conn_4(j)       = connect{1,j}(2,3);
    conn_5(j)       = connect{1,j}(4,3);
end

modulation_1 = fitlm(vector_cov1,modulation_1);
modulation_1 = modulation_1.Residuals.Raw;

modulation_2 = fitlm(vector_cov1,modulation_2);
modulation_2 = modulation_2.Residuals.Raw;

conn_1 = fitlm(vector_cov1,conn_1);
conn_1 = conn_1.Residuals.Raw;

conn_2 = fitlm(vector_cov1,conn_2);
conn_2 = conn_2.Residuals.Raw;

conn_3 = fitlm(vector_cov1,conn_3);
conn_3 = conn_3.Residuals.Raw;

conn_4 = fitlm(vector_cov1,conn_4);
conn_4 = conn_4.Residuals.Raw;

conn_5 = fitlm(vector_cov1,conn_5);
conn_5 = conn_5.Residuals.Raw;

pos_modul_conn_1 = conn_1 + modulation_1;
pos_modul_conn_2 = conn_2 + modulation_2;
neg_modul_conn_1 = conn_1 - modulation_1;
neg_modul_conn_2 = conn_2 - modulation_2;

fn=1;
figure(fn);clf;

red = [1,0.33,0];
grey = [0.5 0.5 0.5];

hold on

% subplot(5, 2, 1);
% v = violinplot(pos_modul_conn_1)
% scatter(1,pos_modul_conn_1(1), 80, [1 0 0], 'filled');
% v.ViolinColor = red;
% v.ShowData = 0;
% [jitter jitterstrength] = parameters_from_violin_to_scatter(pos_modul_conn_1(2:end));
% scatter((1 + jitter.*jitterstrength), pos_modul_conn_1(2:end), 60, grey, 'filled');
% [result(1) p(1)] = ttest2(pos_modul_conn_1(1),pos_modul_conn_1(2:end));
% 
% subplot(5, 2, 2);
% v = violinplot(pos_modul_conn_2)
% scatter(1,pos_modul_conn_2(1), 80, [1 0 0], 'filled');
% v.ViolinColor = red;
% v.ShowData = 0;
% [jitter jitterstrength] = parameters_from_violin_to_scatter(pos_modul_conn_2(2:end));
% scatter((1 + jitter.*jitterstrength), pos_modul_conn_2(2:end), 60, grey, 'filled');
% [result(2) p(2)] = ttest2(pos_modul_conn_2(1),pos_modul_conn_2(2:end));
% 
% 
% subplot(5, 2, 3);
% v = violinplot(neg_modul_conn_1)
% scatter(1,neg_modul_conn_1(1), 80, [1 0 0], 'filled');
% v.ViolinColor = red;
% v.ShowData = 0;
% [jitter jitterstrength] = parameters_from_violin_to_scatter(neg_modul_conn_1(2:end));
% scatter((1 + jitter.*jitterstrength), neg_modul_conn_1(2:end), 60, grey, 'filled');
% [result(3) p(3)] = ttest2(neg_modul_conn_1(1),neg_modul_conn_1(2:end));
% 
% 
% subplot(5, 2, 4);
% v = violinplot(neg_modul_conn_2)
% scatter(1,neg_modul_conn_2(1), 80, [1 0 0], 'filled');
% v.ViolinColor = red;
% v.ShowData = 0;
% [jitter jitterstrength] = parameters_from_violin_to_scatter(neg_modul_conn_2(2:end));
% scatter((1 + jitter.*jitterstrength), neg_modul_conn_2(2:end), 60, grey, 'filled');
% [result(4) p(4)] = ttest2(neg_modul_conn_2(1),neg_modul_conn_2(2:end));


subplot(5, 2, 2);
v = violinplot(modulation_1);
scatter(1,modulation_1(1), 80, [1 0 0], 'filled');
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_1(2:end));
scatter((1 + jitter.*jitterstrength), modulation_1(2:end), 60, grey, 'filled');
[result(1) p(1) ci{1} stat{1}] = ttest2(modulation_1(1),modulation_1(2:end));

sd(1) = std(modulation_1(2:end));
mean_mod(1) = mean(modulation_1(2:end));
nb_sd(1)=(modulation_1(1)-mean(modulation_1(2:end)))/sd(1);

subplot(5, 2, 4);
v = violinplot(modulation_2);
scatter(1,modulation_2(1), 80, [1 0 0], 'filled');
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_2(2:end));
scatter((1 + jitter.*jitterstrength), modulation_2(2:end), 60, grey, 'filled');
[result(2) p(2) ci{2} stat{2}] = ttest2(modulation_2(1),modulation_2(2:end));
sd(2) = std(modulation_2(2:end));
mean_mod(2) = mean(modulation_2(2:end));
nb_sd(2)=(modulation_2(1)-mean(modulation_2(2:end)))/sd(2);

subplot(5, 2, 1);
v = violinplot(conn_1);
scatter(1,conn_1(1), 80, [1 0 0], 'filled');
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_1(2:end));
scatter((1 + jitter.*jitterstrength), conn_1(2:end), 60, grey, 'filled');
[result(3) p(3) ci{3} stat{3}] = ttest2(conn_1(1),conn_1(2:end));
sd(3) = std(conn_1(2:end));
mean_mod(3) = mean(conn_1(2:end));
nb_sd(3)=(conn_1(1)-mean(conn_1(2:end)))/sd(3);

subplot(5, 2, 3);
v = violinplot(conn_2);
scatter(1,conn_2(1), 80, [1 0 0], 'filled');
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_2(2:end));
scatter((1 + jitter.*jitterstrength), conn_2(2:end), 60, grey, 'filled');
[result(4) p(4) ci{4} stat{4}] = ttest2(conn_2(1),conn_2(2:end));
sd(4) = std(conn_2(2:end));
mean_mod(4) = mean(conn_2(2:end));
nb_sd(4)=(conn_2(1)-mean(conn_2(2:end)))/sd(4);

subplot(5, 2, 5);
v = violinplot(conn_3);
scatter(1,conn_3(1), 80, [1 0 0], 'filled');
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_3(2:end));
scatter((1 + jitter.*jitterstrength), conn_3(2:end), 60, grey, 'filled');
[result(5) p(5) ci{5} stat{5}] = ttest2(conn_3(1),conn_3(2:end));
sd(5) = std(conn_3(2:end));
mean_mod(5) = mean(conn_3(2:end));
nb_sd(5)=(conn_3(1)-mean(conn_3(2:end)))/sd(5);

subplot(5, 2, 7);
v = violinplot(conn_4);
scatter(1,conn_4(1), 80, [1 0 0], 'filled');
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_4(2:end));
scatter((1 + jitter.*jitterstrength), conn_4(2:end), 60, grey, 'filled');
[result(6) p(6) ci{6} stat{6}] = ttest2(conn_4(1),conn_4(2:end));
sd(6) = std(conn_4(2:end));
mean_mod(6) = mean(conn_4(2:end));
nb_sd(6)=(conn_4(1)-mean(conn_4(2:end)))/sd(6);

subplot(5, 2, 9);
v = violinplot(conn_5);
scatter(1,conn_5(1), 80, [1 0 0], 'filled');
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_5(2:end));
scatter((1 + jitter.*jitterstrength), conn_5(2:end), 60, grey, 'filled');
[result(7) p(7) ci{7} stat{7}] = ttest2(conn_5(1),conn_5(2:end));
sd(7) = std(conn_5(2:end));
mean_mod(7) = mean(conn_5(2:end));
nb_sd(7)=(conn_5(1)-mean(conn_5(2:end)))/sd(7);

hold off