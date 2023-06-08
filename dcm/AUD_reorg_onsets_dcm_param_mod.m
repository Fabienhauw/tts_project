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

% 3 kinds of models: or 1 = all speech cond vs baseline, or 2 = sentences vs
% scrambled, or 3 = all speech cond vs scrambled.
model_kind = 3; % 1, 2, 3 for previous description of model kinds.

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
    scramble_cond   = 5; % the first 4 conditions are speech vs the 5th = scrambled;
    all_lists_cond  = 2; % the first 2 conditions are list of words or PW;
    clear new_onsets;
    new_onsets(:,1)     = vertcat(onsets{1:nb_cond}); % all the onsets
    new_onsets(:,2)     = vertcat(durations{1:nb_cond}); % all the durations
    if model_kind == 1 % speech cond vs baseline 1 1 1 1 0
        new_onsets2(:,1)    = vertcat(onsets{1:all_speech_cond}); % all the onsets for only words
        new_onsets2(:,2)    = vertcat(durations{1:all_speech_cond}); % all the durations for only words
    elseif model_kind == 2 % sentences vs scrambled 0 0 0 1 -1
        new_onsets2(:,1)    = vertcat(onsets{all_speech_cond}); % all the onsets for only words
        new_onsets2(:,2)    = vertcat(durations{all_speech_cond}); % all the durations for only words
        new_onsets2(:,3)    = vertcat(onsets{scramble_cond}); % all the onsets for only words
        new_onsets2(:,4)    = vertcat(durations{scramble_cond}); % all the durations for only words
    elseif model_kind == 3 % all speech cond vs scrambled 1 1 1 1 -1
        new_onsets2(:,1)    = vertcat(onsets{1:all_speech_cond}); % all the onsets for only words
        new_onsets2(:,2)    = vertcat(durations{1:all_speech_cond}); % all the durations for only words
        new_onsets2(:,3)    = vertcat(onsets{scramble_cond}); % all the onsets for only words
        new_onsets2(:,4)    = vertcat(durations{scramble_cond}); % all the durations for only words
    end
    
    new_onsets3(:,1)    = onsets{1}; % all the onsets for words
    new_onsets3(:,2)    = durations{1}; % all the durations for words
    new_onsets4(:,1)    = onsets{2}; % all the onsets for pw
    new_onsets4(:,2)    = durations{2}; % all the durations for pw
    
    for mp = 1:length(new_onsets(:,1))
        
        if model_kind == 1 % speech cond vs baseline 1 1 1 1 0
            if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
                new_onsets(mp,3)=1; % when speech
            else
                new_onsets(mp,3)=0;
            end
        elseif model_kind == 2 % sentences vs scrambled 0 0 0 1 -1
            if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
                new_onsets(mp,3)=1; % when sent
            elseif ~isempty(find(new_onsets(mp,1) == new_onsets2(:,3)))
                new_onsets(mp,3)=-1; % when scrambled speech
            else
                new_onsets(mp,3)=0;
            end
        elseif model_kind == 3 % all speech cond vs scrambled 1 1 1 1 -1
            if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
                new_onsets(mp,3)=1; % when sent or speech
            else
                new_onsets(mp,3)=-1; % when scrambled speech
            end
        end
        
        
        % for lexical modulation
        if ~isempty(find(new_onsets(mp,1) == new_onsets3(:,1)))
            new_onsets(mp,5)=1; % when words
            new_onsets(mp,6)=-1;
        elseif ~isempty(find(new_onsets(mp,1) == new_onsets4(:,1)))
            new_onsets(mp,5)=-1;
            new_onsets(mp,6)=1;
        else
            new_onsets(mp,5)=0;
            new_onsets(mp,6)=0;
        end
    end
    
    new_onsets  = sortrows(new_onsets);
  
    names2              = 'sounds';
    onsets2{1}          = new_onsets(:,1);
    durations2{1}       = new_onsets(:,2);
    parametric_modul{1} = new_onsets(:,3); % 1 = speech, -1 = scr;
    parametric_modul{2} = new_onsets(:,5); % 1 = words, -1 = pw;

    names=names2;
    onsets=onsets2;
    durations=durations2;
    
    cd(path_to_cpt{k})
    if model_kind == 1 % speech cond vs baseline 1 1 1 1 0
        save('onsets_dcm_param_modul_speech_baseline', 'names', 'onsets', 'durations', 'parametric_modul')
    elseif model_kind == 2 % sentences vs scrambled 0 0 0 1 -1
        save('onsets_dcm_param_modul_sent_scramble', 'names', 'onsets', 'durations', 'parametric_modul')
    elseif model_kind == 3 % all speech cond vs scrambled 1 1 1 1 -1
        save('onsets_dcm_param_modul', 'names', 'onsets', 'durations', 'parametric_modul')
    end
end