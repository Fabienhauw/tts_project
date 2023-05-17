% scripts of a full factorial ANOVA, as seen with Cecile Gallea 2023.
clear;
clc;

res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/Essai_withCecile';

tts_group.path_to_subject;

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet')));
C_droit = S;
C_droit(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
S_droit = S;
S_droit(mask_gauch) = [];

S_effect = [S_droit ; C_droit];

ncon = 5;
levels = 2;
levels_names = {'Group', 'Conditions'};
nlevel(1) = 2;
nlevel(2) = ncon;
count = 0;

for subj = 1 : length(S_droit)
    for con = 1 : ncon
        sub = S_droit(subj).name;
        P{1, con}{subj,1} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,con);
    end
end

for subj = 1 : length(C_droit)
    for con = 1 : ncon
        sub = C_droit(subj).name;
        P{2, con}{subj,1} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,con);
    end
end




matlabbatch{1}.spm.stats.factorial_design.dir = {res_dir};
for tmpLvl = 1 : levels
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).name = levels_names{tmpLvl};
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).levels = nlevel(tmpLvl);
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(tmpLvl).ancova = 0;
end

for tmpLvl1 = 1 : nlevel(1)
    for tmpLvl2 = 1 : nlevel(2)
        count = count + 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(count).levels = [
            tmpLvl1
            tmpLvl2
            ];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(count).scans = P{tmpLvl1,tmpLvl2};
    end
end

matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S_effect(1).name);
exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, S_effect(end).name);
matlabbatch{1}.spm.stats.factorial_design.masking.em = {exp_mask};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(res_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch)

%% design contrasts
words           = [1 0 0 0 0];
pseudowords     = [0 1 0 0 0];
numbers         = [0 0 1 0 0];
normal_speech   = [0 0 0 1 0];
scramble_speech = [0 0 0 0 1];
EOI             = [eye(5)];

phonology       = normal_speech - scramble_speech;
lexicality      = words - pseudowords;

Xvalues = {...
    words, pseudowords, numbers, normal_speech, scramble_speech,...
    lexicality, -lexicality, (words + normal_speech + numbers) - 3*pseudowords,...
    phonology, -phonology, 2*numbers - (words + pseudowords), ...
    normal_speech - words, words + pseudowords + numbers + normal_speech, ...
    words - normal_speech, words - scramble_speech, pseudowords - scramble_speech, (words+pseudowords) - 2*scramble_speech, ...
    (words + normal_speech) - 2*scramble_speech, (normal_speech + pseudowords) - 2*scramble_speech, ...
    (normal_speech + pseudowords + words) - 3*scramble_speech, (normal_speech + pseudowords + words + numbers) - 4*scramble_speech,...
    numbers - words, words - numbers,...
    EOI, ...
    };

Xnames = {...
    'words', 'pseudowords', 'numbers', 'normal_speech', 'scramble_speech',...
    'lexicality', '-lexicality','(words + normal_speech + numbers) - pseudowords',...
    'phonology', '-phonology', 'numbers - (words + pseudowords)', ...
    'normal_speech - words', 'words + pseudowords + numbers + normal_speech',...
    'words - normal_speech', 'words - scramble_speech', 'pseudowords - scramble_speech', '(words+pseudowords) -  scramble_speech', ...
    '(words + normal_speech) -  scramble_speech', '(normal_speech + pseudowords) -  scramble_speech', '(normal_speech + pseudowords + words) -  scramble_speech',...
    '(normal_speech + pseudowords + words + numbers) - scramble_speech',...
    'numbers - words', 'words - numbers', ...
    'EOI', ...
    };

numcomp=length(Xvalues);

% synesthetes
for u=1:numcomp
    values{u}               = [Xvalues{u} zeros(size(Xvalues{u}))];
%     values{u}               = reshape( [repmat(Xvalues{u}, nsub(1),1); repmat(zeros(size(Xvalues{u})), nsub(2),1)],     size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{u}                = ['S ' Xnames{u}];
end
% controls
for u=1:numcomp
    values{numcomp+u}       = [zeros(size(Xvalues{u})) Xvalues{u}];
%     values{numcomp+u}       = reshape([repmat(zeros(size(Xvalues{u})), nsub(1),1); repmat(Xvalues{u}, nsub(2),1)],     size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{numcomp+u}        = ['C ' Xnames{u}];
end
% (synesthetes+controls)
for u=1:numcomp
    values{2*numcomp+u}     = [Xvalues{u} Xvalues{u}];
%     values{2*numcomp+u}     = reshape([repmat(Xvalues{u}, nsub(1),1); repmat(Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{2*numcomp+u}      = ['S+C ' Xnames{u}];
end
% synesthetes-controls
for u=1:numcomp
    values{3*numcomp+u}     = [Xvalues{u} -Xvalues{u}];
%     values{3*numcomp+u}     = reshape([repmat(Xvalues{u}, nsub(1),1); repmat(-Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{3*numcomp+u}      = ['S-C ' Xnames{u}];
end
% controls-synesthetes
for u=1:numcomp
    values{4*numcomp+u}     = [-Xvalues{u} Xvalues{u}];    
%     values{4*numcomp+u}     = reshape([repmat(-Xvalues{u}, nsub(1),1); repmat(Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{4*numcomp+u}      = ['C-S ' Xnames{u}];
end

for n = 1 : numcomp-1
    types{n} = 'T';
    types{numcomp+n} = 'T';
    types{2*numcomp+n} = 'T';
    types{3*numcomp+n} = 'T';
    types{4*numcomp+n} = 'T';
end
for n = numcomp
    types{n}='F';
    types{numcomp+n} = 'F';
    types{2*numcomp+n} = 'F';
    types{3*numcomp+n} = 'F';
    types{4*numcomp+n} = 'F';
end

load SPM;
SPM = rmfield (SPM,'xCon');
% SPM.xCon = struct('name', {}, 'STAT', {}, 'c', {}, 'XO', {}, 'iXO', {}, 'X1o', {}, 'eidf', {}, 'Vcon', {}, 'Vspm', {});
for n=1:max(size(names))
    contrast = spm_FcUtil('Set',names{n}, types{n}, 'c', values{n}', SPM.xX.xKXs); % if error, check "cont" variable line 63, if good number of contrasts (including the resting when appropriate)
    SPM.xCon(n,1) = contrast;
end

save SPM SPM;

spm_contrasts(SPM);