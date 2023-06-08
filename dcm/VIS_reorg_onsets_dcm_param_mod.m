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
    path_to_subj1   = fullfile(D, S(k).name,'Vis/loc/cpt_data');
    path_to_cpt     = [path_to_cpt;path_to_subj1];
end


for k = a : b
    cd(path_to_cpt{k})
    timedata = dir('Timedata*without_resting.mat');
    load(timedata.name)
    nb_cond = length(names)-2; %last 2 conditions are odds and motor responses;
    non_word_cond   = 3; % the first 3 conditions are houses / tools / faces;
    words_cond      = 5; % the 5th conditions is words;
    clear new_onsets;
    new_onsets(:,1)     = vertcat(onsets{1:nb_cond}); % all the onsets
    new_onsets(:,2)     = vertcat(durations{1:nb_cond}); % all the durations
    
    new_onsets2(:,1)    = vertcat(onsets{words_cond}); % all the onsets for only words
    new_onsets2(:,2)    = vertcat(durations{words_cond}); % all the durations for only words
    new_onsets3(:,1)    = vertcat(onsets{1:non_word_cond}); % all the onsets for faces/houses/tools
    new_onsets3(:,2)    = vertcat(durations{1:non_word_cond}); % all the durations for faces/houses/tools

    
    for mp = 1:length(new_onsets(:,1))
        % -1 -1 -1 0 1
        if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
            new_onsets(mp,3)=1; % when words
        elseif ~isempty(find(new_onsets(mp,1) == new_onsets3(:,1)))
            new_onsets(mp,3)=-1;
        else
            new_onsets(mp,3)=0;
        end
    end
    
    new_onsets  = sortrows(new_onsets);
  
    names2              = 'stim';
    onsets2{1}          = new_onsets(:,1);
    durations2{1}       = new_onsets(:,2);
    parametric_modul{1} = new_onsets(:,3); % 1 = words, -1 = non words;

    names=names2;
    onsets=onsets2;
    durations=durations2;
    
    cd(path_to_cpt{k})
    save('onsets_dcm_param_modul_words_non_words', 'names', 'onsets', 'durations', 'parametric_modul')
end