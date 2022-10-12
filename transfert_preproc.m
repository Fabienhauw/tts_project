clear
clc

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';
cd (D);
S = dir(D);
mask = ismember({S.name}, {'.', '..','meinfo.mat'});
S(mask) = [];

a = 41; b = 41;


DF = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/cpt_data_fmri';

path1='/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI/\w*/';
path2='/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';
path_to_subj = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';
subj_avail = dir(path_to_subj);
mask = ismember({subj_avail.name}, {'.', '..'});
subj_avail(mask) = [];
subjects={};
for count = 1:length(subj_avail)
    subjects=[subjects ; subj_avail(count).name];
end

cd (path_to_subj);
if ~exist('a','var')
    disp('Please select the first subject.');
    path = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI/waiting_for_input/';
    while ~isempty(regexp(path,path1)) | ~isempty(regexp(path2,path))
        try path = uigetdir;
            if ~isempty(regexp(path,path1)) %to avoid the selection of a subfolder;
                msg = sprintf('\n%s','Please just select the subject.');
                error(msg);
            elseif ~isempty(regexp(path2,path)) %to avoid the selection of a susfolder;
                msg = sprintf('\n%s','Please just select the subject.');
                error(msg);
            end
        catch MyErr
            disp(MyErr.message)
            subjects
        end
    end
    for k = 1:numel(S)
        if strcmp(path,fullfile(D,S(k).name))==1
            a = k;
        end
    end
else
    a = a;
end

cd (path_to_subj);
if ~exist('b','var')
    disp('Please select the second subject.');
    path = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI/waiting_for_input/';
    while ~isempty(regexp(path,path1)) | ~isempty(regexp(path2,path))
        try path = uigetdir;
            if ~isempty(regexp(path,path1)) %to avoid the selection of a subfolder;
                msg = sprintf('\n%s','Please just select the subject.');
                error(msg);
            elseif ~isempty(regexp(path2,path)) %to avoid the selection of a susfolder;
                msg = sprintf('\n%s','Please just select the subject.');
                error(msg);
            end
        catch MyErr
            disp(MyErr.message)
            subjects
        end
    end
    for k = 1:numel(S)
        if strcmp(path,fullfile(D,S(k).name))==1
            b = k;
        end
    end
else
    b = b;
end

jobs = length(a:b);
% seq_regexp = ['T1w';'Localizer_visual'; 'Localizer_audio';'Resting_state'];
for k = a:b
    fprintf('JOB %d/%d, %s... \n', (k-a+1),jobs, S(k).name);
    directory=fullfile(D, S(k).name);
    cd(directory)
    
    anat = dir('*T1w'); if isempty(anat), anat = dir('*T1*p2'); end; if length(anat)>1, anat = anat(end);, end; anat = anat(1).name;
    loc_vis = dir('*Localizer_visual'); if ~isempty(loc_vis), if length(loc_vis)>1, loc_vis = loc_vis(end); end; loc_vis = loc_vis(1).name;end;
    loc_aud = dir('*Localizer_audio'); if ~isempty(loc_aud), if length(loc_aud)>1, loc_aud = loc_aud(end); end; loc_aud = loc_aud(1).name;end;
    rest_state = dir('*RS'); if ~isempty(rest_state), if length(rest_state)>1, rest_state = rest_state(end); end; rest_state = rest_state(1).name;end;
    loc_col = dir('*Localizer_color'); if ~isempty(loc_col), if length(loc_col)>1, loc_col = loc_col(end); end; loc_col = loc_col(1).name;end;
    loc_unfr_col = dir('*Localizer_unframed_color'); if ~isempty(loc_unfr_col), if length(loc_unfr_col)>1, loc_unfr_col = loc_unfr_col(end); end; loc_unfr_col = loc_unfr_col(1).name;end;
    
    
    %% LOC VISUAL: Get the movement regressors + scan
    try
        cd (fullfile(directory, loc_vis))
        param_1 = dir('v*.json');
        J1 = fullfile(directory,loc_vis,param_1(1).name);
        m_1 = dir('wbet*vtde1.nii');
        mask_1 = fullfile(directory,loc_vis,m_1(1).name);
        ted = (dir('*tedana*')); ted = ted(1).name;
        if exist(ted)
            cd (fullfile(directory,loc_vis,ted));
            mov_reg_1 = dir('multiple_reg*.txt');
            M1 = fullfile(directory,loc_vis,ted,mov_reg_1(1).name);
            simage_1_s4 = dir('s4wts_OC.nii');
            simage_1_s5 = dir('s5wts_OC.nii');
            simage_1_s6 = dir('s6wts_OC.nii');
            wimage_1 = dir('wts_OC.nii');
            image_1 = dir('ts_OC.nii');
            S1_s4 = fullfile(directory,loc_vis,ted,simage_1_s4(1).name);
            S1_s5 = fullfile(directory,loc_vis,ted,simage_1_s5(1).name);
            S1_s6 = fullfile(directory,loc_vis,ted,simage_1_s6(1).name);
            W1 = fullfile(directory,loc_vis,ted,wimage_1(1).name);
            I1 = fullfile(directory,loc_vis,ted,image_1(1).name);
        end
    end
    
    %% LOC AUDIO: Get the movement regressors + scan
    try
        if ~isempty(loc_aud)
            cd (fullfile(directory, loc_aud))
            param_2 = dir('*.json');
            J2=fullfile(directory,loc_aud,param_2(1).name);
            m_2 = dir('wbet*vtde1.nii');
            mask_2 = fullfile(directory,loc_aud,m_2(1).name);
            ted = (dir('*tedana*')); ted = ted(1).name;
            if exist(ted)
                cd (fullfile(directory,loc_aud,ted));
                mov_reg_2 = dir('multiple_reg*.txt');
                M2 = fullfile(directory,loc_aud,ted,mov_reg_2(1).name);
                simage_2_s4 = dir('s4wts_OC.nii');
                simage_2_s5 = dir('s5wts_OC.nii');
                simage_2_s6 = dir('s6wts_OC.nii');
                wimage_2 = dir('wts_OC.nii');
                image_2 = dir('ts_OC.nii');
                S2_s4 = fullfile(directory,loc_aud,ted,simage_2_s4(1).name);
                S2_s5 = fullfile(directory,loc_aud,ted,simage_2_s5(1).name);
                S2_s6 = fullfile(directory,loc_aud,ted,simage_2_s6(1).name);
                W2 = fullfile(directory,loc_aud,ted,wimage_2(1).name);
                I2 = fullfile(directory,loc_aud,ted,image_2(1).name);
            end
        end
    end
    
    %% RESTING STATE: Get the movement regressors + scan
    try
        if ~isempty(rest_state)
            cd (fullfile(directory, rest_state))
            param_3 = dir('v*.json');
            J3=fullfile(directory,rest_state,param_3(1).name);
            m_3 = dir('wbet*vtde1.nii');
            mask_3 = fullfile(directory,rest_state,m_3(1).name);
            ted = (dir('*tedana*')); ted = ted(1).name;
            if exist(ted)
                cd (fullfile(directory,rest_state,ted));
                mov_reg_3 = dir('multiple_reg*.txt');
                M3 = fullfile(directory,rest_state,ted,mov_reg_3(1).name);
                simage_3_s4 = dir('s4wts_OC.nii');
                simage_3_s5 = dir('s5wts_OC.nii');
                simage_3_s6 = dir('s6wts_OC.nii');
                wimage_3 = dir('wts_OC.nii');
                image_3 = dir('ts_OC.nii');
                S3_s4 = fullfile(directory,rest_state,ted,simage_3_s4(1).name);
                S3_s5 = fullfile(directory,rest_state,ted,simage_3_s5(1).name);
                S3_s6 = fullfile(directory,rest_state,ted,simage_3_s6(1).name);
                W3 = fullfile(directory,rest_state,ted,wimage_3(1).name);
                I3 = fullfile(directory,rest_state,ted,image_3(1).name);
            end
        end
    end
    
    %% COLOR LOC: Get the movement regressors + scan
    %     if ~isempty(loc_col)
    %         cd (fullfile(directory, loc_col))
    %         param_4 = dir('v*.json');
    %         J4=fullfile(directory,loc_col,param_4(1).name);
    %         m_4 = dir('wbet*vtde1.nii');
    %         mask_4 = fullfile(directory,loc_col,m_4(1).name);
    %         ted = (dir('*tedana*')); ted = ted(1).name;
    %         if exist(ted)
    %             cd (fullfile(directory,loc_col,ted));
    %             mov_reg_4 = dir('multiple_reg*.txt');
    %             M4 = fullfile(directory,loc_col,ted,mov_reg_4(1).name);
    %             simage_4 = dir('s*wts_OC.nii');
    %             wimage_4 = dir('wts_OC.nii');
    %             image_4 = dir('ts_OC.nii');
    %             S4 = fullfile(directory,loc_col,ted,simage_4(1).name);
    %             W4 = fullfile(directory,loc_col,ted,wimage_4(1).name);
    %             I4 = fullfile(directory,loc_col,ted,image_4(1).name);
    %         end
    %     end
    
    %% UNFRAMED COLOR LOC: Get the movement regressors + scan
    try
        if ~isempty(loc_unfr_col)
            cd (fullfile(directory, loc_unfr_col))
            param_5 = dir('v*.json');
            J5=fullfile(directory,loc_unfr_col,param_5(1).name);
            m_5 = dir('wbet*vtde1.nii');
            mask_5 = fullfile(directory,loc_unfr_col,m_5(1).name);
            ted = (dir('*tedana*')); ted = ted(1).name;
            if exist(ted)
                cd (fullfile(directory,loc_unfr_col,ted));
                mov_reg_5 = dir('multiple_reg*.txt');
                M5 = fullfile(directory,loc_unfr_col,ted,mov_reg_5(1).name);
                simage_5_s4 = dir('s4wts_OC.nii');
                simage_5_s5 = dir('s5wts_OC.nii');
                simage_5_s6 = dir('s6wts_OC.nii');
                wimage_5 = dir('wts_OC.nii');
                image_5 = dir('ts_OC.nii');
                S5_s4 = fullfile(directory,loc_unfr_col,ted,simage_5_s4(1).name);
                S5_s5 = fullfile(directory,loc_unfr_col,ted,simage_5_s5(1).name);
                S5_s6 = fullfile(directory,loc_unfr_col,ted,simage_5_s6(1).name);
                W5 = fullfile(directory,loc_unfr_col,ted,wimage_5(1).name);
                I5 = fullfile(directory,loc_unfr_col,ted,image_5(1).name);
            end
        end
    end
    
    %% Get the warp anat
    cd (fullfile(directory, anat))
    wanat= dir('wmv*.nii');
    S6=fullfile(directory,anat,wanat(1).name);
    wgm = dir('wp1*.nii');
    wwm = dir('wp2*.nii');
    wcsf = dir('wp3*.nii');
    S7=fullfile(directory,anat,wgm(1).name);
    S8=fullfile(directory,anat,wwm(1).name);
    S9=fullfile(directory,anat,wcsf(1).name);
    
    %%%%%%%%%%%%%%%%%%%%
    %% Files movement %%
    %%%%%%%%%%%%%%%%%%%%
    T = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
    
    subject     = strsplit(S(k).name,'_');
    if ~isempty(regexp(subject{end},'Session'))|~isempty(regexp(subject{end},'session'))
        subject = subject(end-1);
    else
        subject = subject(end);
    end
    
    %% CPT data:
    cd (DF)
    data_fold           = dir(sprintf('*%s',subject{1}));
    cd (fullfile(DF,data_fold(1).name));
    aud_onsets          = dir('*Aud*mat');
    vis_onsets          = dir('*Vis.mat');
    col_onsets          = dir(sprintf('*%s_Color.mat',subject{1}));
    if isempty(col_onsets), col_onsets = dir(sprintf('*%s_Col.mat',subject{1})); end
    unfr_col_onsets     = dir('*Unframed_Col.mat');
    if ~isempty(loc_aud)
        for k = 1: length(aud_onsets)
            aud_data_cpt{k} = fullfile(DF,data_fold(1).name,aud_onsets(k).name);
        end
    end
    if ~isempty(loc_vis)
        for k = 1 : length(vis_onsets)
            vis_data_cpt{k} = fullfile(DF,data_fold(1).name,vis_onsets(k).name);
        end
    end
    if ~isempty(loc_col)
        for k = 1: length(col_onsets)
            col_data_cpt{k} = fullfile(DF,data_fold(1).name,col_onsets(k).name);
        end
    end
    
    if ~isempty(loc_unfr_col)
        for k = 1 : length(unfr_col_onsets)
            unfr_col_data_cpt{k} = fullfile(DF,data_fold(1).name,unfr_col_onsets(k).name);
        end
    end
    
    cd (directory)
    if ~exist (fullfile(T,subject{1},'Aud/loc'))
        mkdir(fullfile(T,subject{1},'RS/f'));
        mkdir(fullfile(T,subject{1},'RS/wf'));
        mkdir(fullfile(T,subject{1},'RS/swf'));
        mkdir(fullfile(T,subject{1},'RS/stats'));
        mkdir(fullfile(T,subject{1},'RS/param'));
        mkdir(fullfile(T,subject{1},'Vis/loc/wf'));
        mkdir(fullfile(T,subject{1},'Vis/loc/swf'));
        mkdir(fullfile(T,subject{1},'Vis/loc/f'));
        mkdir(fullfile(T,subject{1},'Vis/loc/stats'));
        mkdir(fullfile(T,subject{1},'Vis/loc/param'));
        mkdir(fullfile(T,subject{1},'Vis/loc/mvpa/balanced_data'));
        mkdir(fullfile(T,subject{1},'Vis/loc/cpt_data'));
        mkdir(fullfile(T,subject{1},'Vis/col/wf'));
        mkdir(fullfile(T,subject{1},'Vis/col/swf'));
        mkdir(fullfile(T,subject{1},'Vis/col/f'));
        mkdir(fullfile(T,subject{1},'Vis/col/stats'));
        mkdir(fullfile(T,subject{1},'Vis/col/param'));
        mkdir(fullfile(T,subject{1},'Vis/col/mvpa/balanced_data'));
        mkdir(fullfile(T,subject{1},'Vis/col/cpt_data'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/wf'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/swf'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/f'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/stats'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/param'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/mvpa/balanced_data'));
        mkdir(fullfile(T,subject{1},'Vis/unfr_col/cpt_data'));
        mkdir(fullfile(T,subject{1},'anat'));
        mkdir(fullfile(T,subject{1},'Aud/loc/wf'));
        mkdir(fullfile(T,subject{1},'Aud/loc/swf'));
        mkdir(fullfile(T,subject{1},'Aud/loc/f'));
        mkdir(fullfile(T,subject{1},'Aud/loc/stats'));
        mkdir(fullfile(T,subject{1},'Aud/loc/param'));
        mkdir(fullfile(T,subject{1},'Aud/loc/mvpa/balanced_data'));
        mkdir(fullfile(T,subject{1},'Aud/loc/cpt_data'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/wf'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/swf'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/f'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/stats'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/param'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/mvpa/balanced_data'));
        mkdir(fullfile(T,subject{1},'Aud/length_task/cpt_data'));
        mkdir(fullfile(T,subject{1},'DTI'));
    end
    
    %% Movement regressors
    %% RESTING STATE
    if ~isempty(rest_state)
        if exist(fullfile(directory,rest_state,ted))==7
            fprintf('copying resting state... ');
            cd (fullfile(T,subject{1},'RS/param'));
            try
                if isempty(dir('*.txt'))
                    copyfile(M3, fullfile(T,subject{1},'RS/param',mov_reg_3(1).name));
                end
                if isempty(dir('wbet*.nii'))
                    copyfile(mask_3, fullfile(T,subject{1},'RS/param',m_3(1).name));
                end
                if isempty(dir('*.json'))
                    copyfile(J3, fullfile(T,subject{1},'RS/param',param_3(1).name));
                end
            end
            cd (fullfile(T,subject{1},'RS/f'));
            if isempty(dir('*.nii'))
                copyfile(I3, fullfile(T,subject{1},'RS/f',image_3(1).name));
            end
            cd (fullfile(T,subject{1},'RS/wf'));
            if isempty(dir('*.nii'))
                copyfile(W3, fullfile(T,subject{1},'RS/wf',wimage_3(1).name));
            end
            cd (fullfile(T,subject{1},'RS/swf'));
            if isempty(dir('s4*.nii'))
                copyfile(S3_s4, fullfile(T,subject{1},'RS/swf',simage_3_s4(1).name));
            end
            if isempty(dir('s5*.nii'))
                copyfile(S3_s5, fullfile(T,subject{1},'RS/swf',simage_3_s5(1).name));
            end
            if isempty(dir('s6*.nii'))
                copyfile(S3_s5, fullfile(T,subject{1},'RS/swf',simage_3_s6(1).name));
            end
        end
    end
    
    %% LOC VISUAL
    if exist(fullfile(directory,loc_vis,ted))==7
        fprintf('copying visual localizer... ');
        cd (fullfile(T,subject{1},'Vis/loc/param'));
        try
            if isempty(dir('*.txt'))
                copyfile(M1, fullfile(T,subject{1},'Vis/loc/param',mov_reg_1(1).name));
            end
            if isempty(dir('wbet*.nii'))
                copyfile(mask_1, fullfile(T,subject{1},'Vis/loc/param',m_1(1).name));
            end
            if isempty(dir('*.json'))
                copyfile(J1, fullfile(T,subject{1},'Vis/loc/param',param_1(1).name));
            end
        end
        cd (fullfile(T,subject{1},'Vis/loc/f'));
        if isempty(dir('*.nii'))
            copyfile(I1, fullfile(T,subject{1},'Vis/loc/f',image_1(1).name));
        end
        cd (fullfile(T,subject{1},'Vis/loc/wf'));
        if isempty(dir('*.nii'))
            copyfile(W1, fullfile(T,subject{1},'Vis/loc/wf',wimage_1(1).name));
        end
        cd (fullfile(T,subject{1},'Vis/loc/swf'));
        if isempty(dir('s4*.nii'))
            copyfile(S1_s4, fullfile(T,subject{1},'Vis/loc/swf',simage_1_s4(1).name));
        end
        if isempty(dir('s5*.nii'))
            copyfile(S1_s5, fullfile(T,subject{1},'Vis/loc/swf',simage_1_s5(1).name));
        end
        if isempty(dir('s6*.nii'))
            copyfile(S1_s6, fullfile(T,subject{1},'Vis/loc/swf',simage_1_s6(1).name));
        end
        cd (fullfile(T,subject{1},'Vis/loc/cpt_data'));
        if isempty(dir('*.mat'))
            for k = 1:length(vis_data_cpt)
                copyfile(vis_data_cpt{k}, fullfile(T,subject{1},'Vis/loc/cpt_data',vis_onsets(k).name));
            end
        end
    end
    
    %% LOC AUDIO
    if ~isempty(loc_aud)
        try
            if exist(fullfile(directory,loc_aud,ted))==7
                fprintf('copying audio localizer... ');
                cd (fullfile(T,subject{1},'Aud/loc/param'));
                
                if isempty(dir('*.txt'))
                    copyfile(M2, fullfile(T,subject{1},'Aud/loc/param',mov_reg_2(1).name));
                end
                if isempty(dir('wbet*.nii'))
                    copyfile(mask_2, fullfile(T,subject{1},'Aud/loc/param',m_2(1).name));
                end
                if isempty(dir('*.json'))
                    copyfile(J2, fullfile(T,subject{1},'Aud/loc/param',param_2(1).name));
                end
                
                cd (fullfile(T,subject{1},'Aud/loc/f'));
                if isempty(dir('*.nii'))
                    copyfile(I2, fullfile(T,subject{1},'Aud/loc/f',image_2(1).name));
                end
                cd (fullfile(T,subject{1},'Aud/loc/wf'));
                if isempty(dir('*.nii'))
                    copyfile(W2, fullfile(T,subject{1},'Aud/loc/wf',wimage_2(1).name));
                end
                cd (fullfile(T,subject{1},'Aud/loc/swf'));
                if isempty(dir('s4*.nii'))
                    copyfile(S2_s4, fullfile(T,subject{1},'Aud/loc/swf',simage_2_s4(1).name));
                end
                if isempty(dir('s5*.nii'))
                    copyfile(S2_s5, fullfile(T,subject{1},'Aud/loc/swf',simage_2_s5(1).name));
                end
                if isempty(dir('s6*.nii'))
                    copyfile(S2_s6, fullfile(T,subject{1},'Aud/loc/swf',simage_2_s6(1).name));
                end
                cd (fullfile(T,subject{1},'Aud/loc/cpt_data'));
                if isempty(dir('*.mat'))
                    for k = 1:length(aud_data_cpt)
                        copyfile(aud_data_cpt{k}, fullfile(T,subject{1},'Aud/loc/cpt_data',aud_onsets(k).name));
                    end
                end
            end
        end
    end
    
    %     if ~isempty(loc_col)
    %         if exist(fullfile(directory,loc_col,ted))==7
    %             cd (fullfile(T,subject{1},'Vis/col/param'));
    %             if isempty(dir('*.txt'))
    %                 copyfile(M4, fullfile(T,subject{1},'Vis/col/param',mov_reg_4(1).name));
    %             end
    %             if isempty(dir('wbet*.nii'))
    %                 copyfile(mask_4, fullfile(T,subject{1},'Vis/col/param',m_4(1).name));
    %             end
    %             if isempty(dir('*.json'))
    %                 copyfile(J4, fullfile(T,subject{1},'Vis/col/param',param_4(1).name));
    %             end
    %             cd (fullfile(T,subject{1},'Vis/col/f'));
    %             if isempty(dir('*.nii'))
    %                 fprintf('copying color localizer... ');
    %                 copyfile(I4, fullfile(T,subject{1},'Vis/col/f',image_4(1).name));
    %             end
    %             cd (fullfile(T,subject{1},'Vis/col/wf'));
    %             if isempty(dir('*.nii'))
    %                 copyfile(W4, fullfile(T,subject{1},'Vis/col/wf',wimage_4(1).name));
    %             end
    %             cd (fullfile(T,subject{1},'Vis/col/swf'));
    %             if isempty(dir('s4*.nii'))
    %                 copyfile(S4, fullfile(T,subject{1},'Vis/col/swf',simage_4(1).name));
    %             end
    %         end
    %         cd (fullfile(T,subject{1},'Vis/col/cpt_data'));
    %         if isempty(dir('*.mat'))
    %             for k = 1:length(col_data_cpt)
    %                 copyfile(col_data_cpt{k}, fullfile(T,subject{1},'Vis/col/cpt_data',col_onsets(k).name));
    %             end
    %         end
    %     end
    
    %% LOC UNFR COLOR
    if ~isempty(loc_unfr_col)
        if exist(fullfile(directory,loc_unfr_col,ted))==7
            fprintf('copying unframed color localizer... \n');
            cd (fullfile(T,subject{1},'Vis/unfr_col/param'));
            try
                if isempty(dir('*.txt'))
                    copyfile(M5, fullfile(T,subject{1},'Vis/unfr_col/param',mov_reg_5(1).name));
                end
                if isempty(dir('wbet*.nii'))
                    copyfile(mask_5, fullfile(T,subject{1},'Vis/unfr_col/param',m_5(1).name));
                end
                if isempty(dir('*.json'))
                    copyfile(J5, fullfile(T,subject{1},'Vis/unfr_col/param',param_5(1).name));
                end
            end
            cd (fullfile(T,subject{1},'Vis/unfr_col/f'));
            if isempty(dir('*.nii'))
                copyfile(I5, fullfile(T,subject{1},'Vis/unfr_col/f',image_5(1).name));
            end
            cd (fullfile(T,subject{1},'Vis/unfr_col/wf'));
            if isempty(dir('*.nii'))
                copyfile(W5, fullfile(T,subject{1},'Vis/unfr_col/wf',wimage_5(1).name));
            end
            cd (fullfile(T,subject{1},'Vis/unfr_col/swf'));
            if isempty(dir('s4*.nii'))
                copyfile(S5_s4, fullfile(T,subject{1},'Vis/unfr_col/swf',simage_5_s4(1).name));
            end
            if isempty(dir('s5*.nii'))
                copyfile(S5_s5, fullfile(T,subject{1},'Vis/unfr_col/swf',simage_5_s5(1).name));
            end
            if isempty(dir('s6*.nii'))
                copyfile(S5_s6, fullfile(T,subject{1},'Vis/unfr_col/swf',simage_5_s6(1).name));
            end
        end
        cd (fullfile(T,subject{1},'Vis/unfr_col/cpt_data'));
        if isempty(dir('*.mat'))
            for k = 1:length(unfr_col_data_cpt)
                copyfile(unfr_col_data_cpt{k}, fullfile(T,subject{1},'Vis/unfr_col/cpt_data',unfr_col_onsets(k).name));
            end
        end
    end
    
    cd (fullfile(T,subject{1},'anat'));
    if ~exist(wanat(1).name, 'file')
        copyfile (S6, fullfile(T,subject{1},'anat'));
        copyfile (S7, fullfile(T,subject{1},'anat'));
        copyfile (S8, fullfile(T,subject{1},'anat'));
        copyfile (S9, fullfile(T,subject{1},'anat'));
    end
end

fprintf('DONE! \n');

% cd '/home/fabien/Desktop/Fabien/NeoTopLex'
% Batches.choose_subj.stats.timedata_reorg;
%
% cd '/home/fabien/Desktop/Fabien/NeoTopLex'
% Batches.choose_subj.stats.import_timedata;
%
% cd '/home/fabien/Desktop/Fabien/NeoTopLex'