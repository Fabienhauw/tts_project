%%
clear;clc;

addpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12')
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/MATLAB/spm12/matlabbatch'))

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

a = 1; b = 48;

path_to_cpt = {};

for k = a : b
    path_to_subj1   = fullfile(D, S(k).name,'Aud/loc/cpt_data');
    path_to_cpt     = [path_to_cpt;path_to_subj1];
end


for k = a : b
    cd(path_to_cpt{k})
    timedata = dir('Timedata*without_resting.mat');
    load(timedata.name)
    nb_cond = length(names)-2; %last 2 conditions are odds and motor responses, to check with Jean D.;
    all_speech_cond = 4; % the first 4 conditions are speech vs the 5th = scrambled;
    all_lists_cond  = 2; % the first 2 conditions are list of words or PW;
    clear new_onsets;
    new_onsets(:,1)     = vertcat(onsets{1:nb_cond}); % all the onsets
    new_onsets(:,2)     = vertcat(durations{1:nb_cond}); % all the durations
    new_onsets2(:,1)    = vertcat(onsets{1:all_speech_cond}); % all the onsets for only words
    new_onsets2(:,2)    = vertcat(durations{1:all_speech_cond}); % all the onsets for only words
    
    new_onsets3(:,1)    = vertcat(onsets{1:all_lists_cond}); % all the onsets for words and PW
    new_onsets3(:,2)    = vertcat(durations{1:all_lists_cond}); % all the durations for words and PW
    new_onsets4(:,1)    = onsets{1}; % all the onsets for words
    new_onsets4(:,2)    = durations{1}; % all the durations for words
    
    for mp = 1:length(new_onsets(:,1))
        if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
            new_onsets(mp,3)=1;
            new_onsets(mp,4)=-1;
        else
            new_onsets(mp,3)=-1;
            new_onsets(mp,4)=1;
        end
    end
    
    for mp = 1:length(new_onsets3(:,1))
        if ~isempty(find(new_onsets3(mp,1) == new_onsets4(:,1)))
            new_onsets3(mp,3)=1;
            new_onsets3(mp,4)=-1;
        else
            new_onsets3(mp,3)=-1;
            new_onsets3(mp,4)=1;
        end
    end
    
    new_onsets  = sortrows(new_onsets);
    new_onsets3 = sortrows(new_onsets3);
  
    names2              = {'sounds'; 'lists'};
    onsets2{1}          = new_onsets(:,1);
    onsets2{2}          = new_onsets3(:,1);
    durations2{1}       = new_onsets(:,2);
    durations2{2}       = new_onsets3(:,2);
    parametric_modul{1} = new_onsets(:,3);
    parametric_modul{2} = new_onsets3(:,3);

    names=names2;
    onsets=onsets2;
    durations=durations2;
    
    cd(path_to_cpt{k})
    save('onsets_dcm_param_modul', 'names', 'onsets', 'durations', 'parametric_modul')
end