clear;clc;
addpath('/home/fabien.hauw/MATLAB Add-Ons/spm12');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers/DCM');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers');
addpath('/home/fabien.hauw/Desktop/Fabien/NeoTopLex/Divers/plot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   uncomment to apply the scripts to all subjects/controls   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;
Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..','FLC_func'});
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

% get coordinates for all VOI used in DCM for all subjects

Batches.path_to_subject;
mask = ismember({S.name}, {'.', '..', 'FLC_func'});
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
for k = 1:nb_subj
    if exist(fullfile(path_to_scan{k}, 'Aud/stats/chaperon/stc'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/stc/dcm_model_param_modul_endroit');
    else
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit');
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

new_nsubj = nb_subj - sum(mask);
path_to_scan(mask) = '';

for k=1:new_nsubj
    if isempty (regexp(path_to_scan{k}, 'FLC'))
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit/dcm/all_models');
    else
        path_to_stats=fullfile(path_to_scan{k}, 'Aud/stats/run/dcm_model_param_modul_endroit/dcm_minus_Sujet08/all_models');
    end
    cd(path_to_stats)
    dcm = dir('DCM*struct_9_num_56*');
    load(dcm.name)

    connect{k} = DCM.Ep.A;
    modul{k} = DCM.Ep.B;
end

modulation_1 = [];
modulation_2 = [];
modulation_3 = [];
modulation_4 = [];
modulation_5 = [];

conn_1 = [];
conn_2 = [];
conn_3 = [];
conn_4 = [];
conn_5 = [];
conn_6 = [];

for j=1:new_nsubj
    modulation_1(j) = modul{1,j}(2,1,2);
    modulation_2(j) = modul{1,j}(3,1,2);
    modulation_3(j) = modul{1,j}(4,2,2);
    modulation_4(j) = modul{1,j}(2,3,2);
    modulation_5(j) = modul{1,j}(4,3,2);
    
    conn_1(j)       = connect{1,j}(2,1);
    conn_2(j)       = connect{1,j}(3,1);
    conn_3(j)       = connect{1,j}(3,2);
    conn_4(j)       = connect{1,j}(4,2);
    conn_5(j)       = connect{1,j}(2,3);
    conn_6(j)       = connect{1,j}(4,3);
end

pos_modul_conn_1 = conn_1 + modulation_1;
pos_modul_conn_2 = conn_2 + modulation_2;
pos_modul_conn_3 = conn_4 + modulation_3;
pos_modul_conn_4 = conn_5 + modulation_4;
pos_modul_conn_5 = conn_6 + modulation_5;
neg_modul_conn_1 = conn_1 - modulation_1;
neg_modul_conn_2= conn_2 - modulation_2;
neg_modul_conn_3= conn_4 - modulation_3;
neg_modul_conn_4= conn_5- modulation_4;
neg_modul_conn_5 = conn_6- modulation_5;

fn=1;
figure(fn);clf;
hold on

% subplot(5, 2, 1);
% violinplot(pos_modul_conn_1)
% scatter(1,pos_modul_conn_1(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 2);
% violinplot(pos_modul_conn_2)
% scatter(1,pos_modul_conn_2(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 3);
% violinplot(pos_modul_conn_3)
% scatter(1,pos_modul_conn_3(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 4);
% violinplot(pos_modul_conn_4)
% scatter(1,pos_modul_conn_4(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 5);
% violinplot(neg_modul_conn_1)
% scatter(1,neg_modul_conn_1(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 6);
% violinplot(neg_modul_conn_2)
% scatter(1,neg_modul_conn_2(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 7);
% violinplot(neg_modul_conn_3)
% scatter(1,neg_modul_conn_3(1), 80, [1 0 0], 'filled')
% 
% subplot(5, 2, 8);
% violinplot(neg_modul_conn_4)
% scatter(1,neg_modul_conn_4(1), 80, [1 0 0], 'filled')

red = [1,0.33,0];
grey = [0.5 0.5 0.5];

modulation_1    = modulation_1';
modulation_2    = modulation_2';
modulation_3    = modulation_3';
modulation_4    = modulation_4';
modulation_5    = modulation_5';
modulation      = [modulation_1 modulation_2 modulation_3 modulation_4 modulation_5];

i=1;
subplot(6, 2, 2);
v = violinplot(modulation_1);
scatter(1,modulation_1(1), 80, [1 0 0], 'filled')
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_1(2:end));
scatter((1 + jitter.*jitterstrength), modulation_1(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(modulation_1(1),modulation_1(2:end));
sd(i) = std(modulation_1(2:end));
mean_mod(i) = mean(modulation_1(2:end));
nb_sd(i)=(modulation_1(1)-mean(modulation_1(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 4);
v = violinplot(modulation_2);
scatter(1,modulation_2(1), 80, [1 0 0], 'filled')
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_2(2:end));
scatter((1 + jitter.*jitterstrength), modulation_2(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(modulation_2(1),modulation_2(2:end));
sd(i) = std(modulation_2(2:end));
mean_mod(i) = mean(modulation_2(2:end));
nb_sd(i)=(modulation_2(1)-mean(modulation_2(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 8);
v = violinplot(modulation_3);
scatter(1,modulation_3(1), 80, [1 0 0], 'filled')
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_3(2:end));
scatter((1 + jitter.*jitterstrength), modulation_3(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(modulation_3(1),modulation_3(2:end));
sd(i) = std(modulation_3(2:end));
mean_mod(i) = mean(modulation_3(2:end));
nb_sd(i)=(modulation_3(1)-mean(modulation_3(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 10);
v = violinplot(modulation_4);
scatter(1,modulation_4(1), 80, [1 0 0], 'filled')
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_4(2:end));
scatter((1 + jitter.*jitterstrength), modulation_4(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(modulation_4(1),modulation_4(2:end));
sd(i) = std(modulation_4(2:end));
mean_mod(i) = mean(modulation_4(2:end));
nb_sd(i)=(modulation_4(1)-mean(modulation_4(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 12);
v = violinplot(modulation_5);
scatter(1,modulation_5(1), 80, [1 0 0], 'filled')
v.ViolinColor = red;
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(modulation_5(2:end));
scatter((1 + jitter.*jitterstrength), modulation_5(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(modulation_5(1),modulation_5(2:end));
sd(i) = std(modulation_5(2:end));
mean_mod(i) = mean(modulation_5(2:end));
nb_sd(i)=(modulation_5(1)-mean(modulation_5(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 1);
v = violinplot(conn_1);
scatter(1,conn_1(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_1(2:end));
scatter((1 + jitter.*jitterstrength), conn_1(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_1(1),conn_1(2:end));
sd(i) = std(conn_1(2:end));
mean_mod(i) = mean(conn_1(2:end));
nb_sd(i)=(conn_1(1)-mean(conn_1(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 3);
v = violinplot(conn_2);
scatter(1,conn_2(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_2(2:end));
scatter((1 + jitter.*jitterstrength), conn_2(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_2(1),conn_2(2:end));
sd(i) = std(conn_2(2:end));
mean_mod(i) = mean(conn_2(2:end));
nb_sd(i)=(conn_2(1)-mean(conn_2(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 5);
v = violinplot(conn_3);
scatter(1,conn_3(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_3(2:end));
scatter((1 + jitter.*jitterstrength), conn_3(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_3(1),conn_3(2:end));
sd(i) = std(conn_3(2:end));
mean_mod(i) = mean(conn_3(2:end));
nb_sd(i)=(conn_3(1)-mean(conn_3(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 7);
v = violinplot(conn_4);
scatter(1,conn_4(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_4(2:end));
scatter((1 + jitter.*jitterstrength), conn_4(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_4(1),conn_4(2:end));
sd(i) = std(conn_4(2:end));
mean_mod(i) = mean(conn_4(2:end));
nb_sd(i)=(conn_4(1)-mean(conn_4(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 9);
v = violinplot(conn_5);
scatter(1,conn_5(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_5(2:end));
scatter((1 + jitter.*jitterstrength), conn_5(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_5(1),conn_5(2:end));
sd(i) = std(conn_5(2:end));
mean_mod(i) = mean(conn_5(2:end));
nb_sd(i)=(conn_5(1)-mean(conn_5(2:end)))/sd(i);

i = i + 1;
subplot(6, 2, 11);
v = violinplot(conn_6);
scatter(1,conn_6(1), 80, [1 0 0], 'filled')
v.ShowData = 0;
[jitter jitterstrength] = parameters_from_violin_to_scatter(conn_6(2:end));
scatter((1 + jitter.*jitterstrength), conn_6(2:end), 60, grey, 'filled');
[result(i) p(i) ci{i} stat{i}] = ttest2(conn_6(1),conn_6(2:end));
sd(i) = std(conn_6(2:end));
mean_mod(i) = mean(conn_6(2:end));
nb_sd(i)=(conn_6(1)-mean(conn_6(2:end)))/sd(i);

hold off