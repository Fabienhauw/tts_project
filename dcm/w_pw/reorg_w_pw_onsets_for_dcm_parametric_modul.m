%%
Batches.path_to_subject;
S2=S;
mask2 = ismember({S.name}, {'FLC_func'});
S2(mask2) = [];
subjects1={S2.name};
path_to_stats = {};
path_to_cpt = {};
path_to_json = {};
for k=1:numel(S)
    path_to_subj1 = fullfile(D, subjects1{k},'Aud/cpt_data/run');
    path_to_subj2 = fullfile(D, subjects1{k},'Aud/stats/run');
    path_to_subj3 = fullfile(D, subjects1{k},'Aud/nb_vol/run');
    path_to_cpt = [path_to_cpt;path_to_subj1];
    path_to_stats = [path_to_stats;path_to_subj2];
    path_to_json = [path_to_json;path_to_subj3];
end

subjects=[subjects1];

for k=1:length(subjects)
    cd(path_to_cpt{k})
    timedata=dir('Timedata*Aud.mat');
    load(timedata.name)
    nb_cond=length(names)-2; %-2 because the last 2 conditions are odds and motor responses, only words and pw interest us in dcm;
    allwords_cond=4; % the 4 firsts conditions are words;
    clear new_onsets;
    new_onsets(:,1)=vertcat(onsets{1:nb_cond}); %all the onsets
    new_onsets(:,2)=vertcat(durations{1:nb_cond}); %all the durations
    new_onsets2(:,1)=vertcat(onsets{1:allwords_cond}); %all the onsets for only words
    new_onsets2(:,2)=vertcat(durations{1:allwords_cond}); %all the onsets for only words
    
    for mp = 1:length(new_onsets(:,1))
        if ~isempty(find(new_onsets(mp,1) == new_onsets2(:,1)))
            new_onsets(mp,3)=1;
            new_onsets(mp,4)=-1;
        else
            new_onsets(mp,3)=-1;
            new_onsets(mp,4)=1;
        end
    end
    
    new_onsets  = sortrows(new_onsets);
    new_onsets2 = sortrows(new_onsets2); % sort the matrice used for words only;
        
    cd (path_to_json{k});
    json=dir('*.json');
    json=json.name;
    res = get_string_from_json(json, {'CsaSeries.MrPhoenixProtocol.lRepetitions', 'RepetitionTime'}, {'num', 'num'});
    nbVol       = res{1}+1;
    TR          = res{2}/1000; 
    total_time=nbVol*TR;
        
    names2{1}='sounds';
    onsets2{1}          = new_onsets(:,1);
    durations2{1}       = new_onsets(:,2);
    parametric_modul    = new_onsets(:,3);
    parametric_modul_2  = new_onsets(:,4);

    names=names2;
    onsets=onsets2;
    durations=durations2;
    
    cd(path_to_cpt{k})
    save('onsets_run_dcm_param_modul', 'names', 'onsets', 'durations', 'parametric_modul', 'parametric_modul_2')
end