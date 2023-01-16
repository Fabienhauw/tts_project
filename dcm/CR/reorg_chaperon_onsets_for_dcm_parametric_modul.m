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
    path_to_subj1 = fullfile(D, subjects1{k},'Aud/cpt_data/chaperon');
    path_to_subj2 = fullfile(D, subjects1{k},'Aud/stats/chaperon');
    path_to_subj3 = fullfile(D, subjects1{k},'Aud/nb_vol/chaperon');
    path_to_cpt = [path_to_cpt;path_to_subj1];
    path_to_stats = [path_to_stats;path_to_subj2];
    path_to_json = [path_to_json;path_to_subj3];
end

Batches.rosso.path_to_subject;
S3=S;
subjects2={S3.name};
for k=1:numel(S)
    path_to_subj1 = fullfile(D, subjects2{k},'Aud/cpt_data/chaperon');
    path_to_subj2 = fullfile(D, subjects2{k}, 'Aud/stats/chaperon/stc');
    path_to_subj3 = fullfile(D, subjects2{k},'Aud/nb_vol/chaperon');
    path_to_cpt = [path_to_cpt;path_to_subj1];
    path_to_stats = [path_to_stats;path_to_subj2];
    path_to_json = [path_to_json;path_to_subj3];
end

subjects=[subjects1 subjects2];

for k=1:length(subjects)
    cd(path_to_cpt{k})
    load('onsets_chaperon.mat')
    nb_cond=length(names);
    clear new_onsets;
    new_onsets(:,1)=[onsets{1:nb_cond}]; %all the onsets
    new_onsets(:,2)=[durations{1:nb_cond}]; %all the durations
    
    for mp = 1:length(new_onsets(:,1))
        if ~isempty(find(new_onsets(mp,1) == onsets{1,1}))
            new_onsets(mp,3)=1;
            new_onsets(mp,4)=-1;
        else
            new_onsets(mp,3)=-1;
            new_onsets(mp,4)=1;
        end
    end
    
    new_onsets  = sortrows(new_onsets);
    
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
    save('onsets_chap_dcm_param_modul', 'names', 'onsets', 'durations', 'parametric_modul','parametric_modul_2')
end