dir1 = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/syn_group_rs_rh_s5_brainnetome/results/preprocessing';
dir2 = '/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/GraphVar_2.03a/workspaces/syn_con_brainnetome/data/Signals';
cd (dir1)
if ~isdir(dir2)
    mkdir(dir2)
end
subjs = dir ('ROI_Subject*Condition001.mat');
for ss = 1 : length(subjs)
    clear ROISignals
    subj_name = regexp(subjs(ss).name, 'Subject\d\d\d', 'match');
    load(subjs(ss).name)
    data(:,250:251) = '';
    data(:,1:3) = '';
    for dd = 1 : size(data,2)
        ROISignals(:,dd) = data{1,dd};
    end
    save(sprintf('%s/ROISignals_%s.mat', dir2, subj_name{1}), 'ROISignals')
end

%%
dir1 = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/rs/syn_group_rs_rh_s5_brainnetome/results/firstlevel/SBC_01';
dir2 = '/home/fabien.hauw/Desktop/toolbox_matlab/GraphVar_2.03a/workspaces/syn_con_rh_brainnetome/data/CorrMatrix';
if ~isdir(dir2)
    mkdir(dir2)
end
cd (dir1)
subjs = dir ('resultsROI_Subject*.mat');
for ss = 1 : length(subjs)
    clear CorrMatrix
    subj_name = regexp(subjs(ss).name, 'Subject\d\d\d', 'match');
    load(subjs(ss).name)
    Z(:,247) = '';
    CorrMatrix = Z;
    save(sprintf('%s/CorrMatrix_%s.mat', dir2, subj_name{1}), 'CorrMatrix')
end
