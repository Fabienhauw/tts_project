clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
cd(D)
S1 = dir('*Control*'); S2 = dir('*Sujet*');
S = [S1;S2];

for k = 1:2
    
    cd (fullfile(D, S(k).name))
    timedata_raw = dir('Timedata*Aud.mat');
    results_raw = dir('result*Aud*');
    
    load(timedata_raw(1).name);
    load(results_raw(1).name);
    
    names(8)            = {'Resting'};
    j                   = 1;
    onsets{8,1}(j,1)      = 0;
    durations{8,1}(j,1)   = pict2(1).start;
    for pp = 1 : length(pict2)
        
        if isequal(pict2(pp).catcontent, 'Resting') & pict2(pp).duration > 1
            j = j+1;
            onsets{8,1}(j,1)        = pict2(pp).start;
            durations{8,1}(j,1)     = pict2(pp).duration;
        end
    end
    
    save(timedata_raw.name,'names','onsets','durations');

end

%% This part is to delete the resting regressor

clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
cd(D)
S1 = dir('*Control*'); S2 = dir('*Sujet*');
S = [S1;S2];

for k = 1 : length(S)
    
    cd (fullfile(D, S(k).name))
    timedata_raw = dir('Timedata*Aud.mat');
    try delete(timedata_raw2.name), catch, end;
    load(timedata_raw(1).name);
    
    names(8)        = ''; %because the last regressor recorded during the experiment was the motor category, and the subject was not asked to click on the button because blind.
    for item = 1 : length(onsets) - 1
        onsets2{item, 1}          = onsets{item};
        durations2{item, 1}       = durations{item};
    end
    clear onsets; clear durations;
    onsets = onsets2; durations = durations2;
    
    timedata_raw.name = [timedata_raw.name(1:end-4) '_without_resting.mat'];
    save(timedata_raw.name,'names','onsets','durations');

end

%% same for visual exp
% clear;clc;
% 
% D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
% cd(D)
% S1 = dir('*Control*'); S2 = dir('*Sujet*');
% S = [S1;S2];
% 
% for k = 1:length(S)
%     
%     cd (fullfile(D, S(k).name))
%     timedata_raw = dir('Timedata*Vis.mat');
%     results_raw = dir('result*Vis*');
%     
%     load(timedata_raw(1).name);
%     load(results_raw(1).name);
%     
%     names(8)            = {'Resting'};
%     j                   = 1;
%     onsets{8,1}(j,1)      = 0;
%     durations{8,1}(j,1)   = pict(1).start;
%     for pp = 1 : length(pict)
%         
%         if isequal(pict(pp).catcontent, 'Resting') & pict(pp).duration > 1
%             j = j+1;
%             onsets{8,1}(j,1)        = pict(pp).start;
%             durations{8,1}(j,1)     = pict(pp).duration;
%         end
%     end
%     
%     save(timedata_raw.name,'names','onsets','durations');
% 
% end
% 
% %% This part is to delete the resting regressor
% 
% clear;clc;
% 
% D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
% cd(D)
% S1 = dir('*Control*'); S2 = dir('*Sujet*');
% S = [S1;S2];
% 
% for k = 1 : length(S)
%     
%     cd (fullfile(D, S(k).name))
%     timedata_raw = dir('Timedata*Vis.mat');
%     try delete(timedata_raw2.name), catch, end;
%     load(timedata_raw(1).name);
%     
%     names(8)        = ''; %because the last regressor recorded during the experiment was the motor category, and the subject was not asked to click on the button because blind.
%     for item = 1 : length(onsets) - 1
%         onsets2{item, 1}          = onsets{item};
%         durations2{item, 1}       = durations{item};
%     end
%     clear onsets; clear durations;
%     onsets = onsets2; durations = durations2;
%     
%     timedata_raw.name = [timedata_raw.name(1:end-4) '_without_resting.mat'];
%     save(timedata_raw.name,'names','onsets','durations');
% 
% end

%% same for color exp
clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
cd(D)
S1 = dir('*Control*'); S2 = dir('*Sujet*');
S = [S1;S2];

for k = 27
    
    cd (fullfile(D, S(k).name))
    timedata_raw = dir('Timedata*Col*'); if length(timedata_raw) >1, timedata_raw = dir('Timedata*Unframed_Col.mat'); end
    results_raw = dir('result*Col*'); if length(results_raw) >1, results_raw = dir('result*Unframed_Col.mat'); end
    
    load(timedata_raw(1).name);
    load(results_raw(1).name);
    
    names(3)            = {'Resting'};
    j                   = 1;
    onsets{3,1}(j,1)      = 0;
    durations{3,1}(j,1)   = pict(1).start;
    for pp = 1 : length(pict)
        
        if isequal(pict(pp).catcontent, 'Resting') & pict(pp).duration > 1
            j = j+1;
            onsets{3,1}(j,1)        = pict(pp).start;
            durations{3,1}(j,1)     = pict(pp).duration;
        end
    end
    
    save(timedata_raw.name,'names','onsets','durations');

end

%% This part is to delete the resting regressor

clear;clc;

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';
cd(D)
S1 = dir('*Control*'); S2 = dir('*Sujet*');
S = [S1;S2];

for k = 1 : length(S)
    
    cd (fullfile(D, S(k).name))
    timedata_raw = dir('Timedata*Col*'); 
    if length(timedata_raw) >1 & length(dir('Timedata*Unframed_Col.mat'))>0
        timedata_raw = dir('Timedata*Unframed_Col.mat');
    elseif length(timedata_raw) >1 & any(endsWith({timedata_raw.name},'without_resting.mat'))
        timedata_raw(endsWith({timedata_raw.name},'without_resting.mat')) = [];
    end
    try delete(timedata_raw2.name), catch, end;
    load(timedata_raw(1).name);
    
    names(3)        = ''; %because the last regressor recorded during the experiment was the motor category, and the subject was not asked to click on the button because blind.
    for item = 1 : length(onsets) - 1
        onsets2{item, 1}          = onsets{item};
        durations2{item, 1}       = durations{item};
    end
    clear onsets; clear durations;
    onsets = onsets2; durations = durations2;
    
    timedata_raw.name = [timedata_raw.name(1:end-4) '_without_resting.mat'];
    save(timedata_raw.name,'names','onsets','durations');

end