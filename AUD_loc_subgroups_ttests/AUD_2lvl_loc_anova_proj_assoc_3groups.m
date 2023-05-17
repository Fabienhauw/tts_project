clear;
clc;

res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ANOVA_s5_proj_assoc_without_resting';

D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
mask = ismember({S.name}, {'.', '..'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con = S;
S_con(mask_gauch_con) = [];
con_group = repmat({'control'}, 1, length(S_con));

mask_proj =  cellfun(@isempty,(regexp({S.name},'Sujet01|Sujet03|Sujet06|Sujet15')));
Syn_proj = S;
Syn_proj(mask_proj) = [];
syn_proj_group = repmat({'syn_proj'}, 1, length(Syn_proj));

mask_asso =  ~cellfun(@isempty,(regexp({S.name},'Sujet03|Sujet01|Sujet15|Sujet06|Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
Syn_asso = S;
Syn_asso(mask_asso) = [];
syn_asso_group = repmat({'syn_asso'}, 1, length(Syn_asso));

S_effect = [Syn_proj ; Syn_asso; S_con];

categ = [syn_proj_group, syn_asso_group, con_group];
nsub(1) = length(Syn_proj); % synesthetes projectors
nsub(2) = length(Syn_asso); % synesthetes associators
nsub(3) = length(S_con); % controls

levels = 3;
levels_names = {'Group', 'Conditions', 'Sujet'};

ngroup = 3;
ncon = 5;
totsub = numel(S_effect);
nscan = ncon*totsub;

nlevel(1) = 2;
nlevel(2) = ncon;
nlevel(3) = totsub;

subname = {S_effect.name};

totsub=length(subname);

%-----------------------------------------------------------------
for subj = 1 : length(Syn_proj)
    for con = 1 : ncon
        sub = Syn_proj(subj).name;
        P{1, subj}{con,1} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,con);
    end
end

for subj = 1 : length(Syn_asso)
    for con = 1 : ncon
        sub = Syn_asso(subj).name;
        P{2, subj}{con,1} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,con);
    end
end

for subj = 1 : length(S_con)
    for con = 1 : ncon
        sub = S_con(subj).name;
        P{3, subj}{con,1} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,con);
    end
end

%-----------------------------------------------------------------

matlabbatch{1}.spm.stats.factorial_design.dir = {res_dir};

% group
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Group';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0; % 0 = independance, 1 = dependance;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1; % 0 = equal, 1 = unequal
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;

% condition
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Conditions';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

% sujet
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Sujet';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0; 
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

count = 0;

for k = 1 : numel(Syn_proj)
    count = count + 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).scans = P{1,k};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).conds = [
        repmat(1, 1, ncon);1:ncon; repmat(count, 1, ncon)]';
end
for k = 1 : numel(Syn_asso)
    count = count + 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).scans = P{2,k};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).conds = [
        repmat(2, 1, ncon);1:ncon; repmat(count, 1, ncon)]';
end
for k = 1 : numel(S_con)
    count = count + 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).scans = P{3,k};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).conds = [
        repmat(3, 1, ncon);1:ncon; repmat(count, 1, ncon)]';
end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [
    1
    2
    ];
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 3;

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

%%
% Contrasts----------------------------------------------------------------

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
n1 = nsub(1); n2 = nsub(2); n3 = nsub(3);

% synesthetes proj
for u=1:numcomp
    values{u}               = [ones(size(Xvalues{u},1),n1)/n1 zeros(size(Xvalues{u},1),n2) zeros(size(Xvalues{u},1),n3) Xvalues{u} zeros(size(Xvalues{u})) zeros(size(Xvalues{u}))];
    names{u}                = ['Sp ' Xnames{u}];
end
% synesthetes assoc
for u=1:numcomp
    values{numcomp+u}       = [zeros(size(Xvalues{u},1),n1) ones(size(Xvalues{u},1),n2)/n2 zeros(size(Xvalues{u},1),n3) zeros(size(Xvalues{u})) Xvalues{u} zeros(size(Xvalues{u}))];
    names{numcomp+u}        = ['Sa ' Xnames{u}];
end

% controls
for u=1:numcomp
    values{2*numcomp+u}       = [zeros(size(Xvalues{u},1),n1) zeros(size(Xvalues{u},1),n2) ones(size(Xvalues{u},1),n3)/n3 zeros(size(Xvalues{u})) zeros(size(Xvalues{u})) Xvalues{u}];
    names{2*numcomp+u}        = ['C ' Xnames{u}];
end

% (synesthetes+controls)
for u=1:numcomp
    values{3*numcomp+u}     = [ones(size(Xvalues{u},1),n1)/n1 ones(size(Xvalues{u},1),n2)/n2 ones(size(Xvalues{u},1),n3)/n3 Xvalues{u} Xvalues{u} Xvalues{u}];
    names{3*numcomp+u}      = ['S+C ' Xnames{u}];
end

% proj - assoc
for u=1:numcomp
    values{4*numcomp+u}     = [ones(size(Xvalues{u},1),n1)/n1 -ones(size(Xvalues{u},1),n2)/n2 zeros(size(Xvalues{u},1),n3) Xvalues{u} -Xvalues{u} zeros(size(Xvalues{u}))];
    names{4*numcomp+u}      = ['Sp-Sa ' Xnames{u}];
end

% assoc -  proj
for u=1:numcomp
    values{5*numcomp+u}     = [-ones(size(Xvalues{u},1),n1)/n1 ones(size(Xvalues{u},1),n2)/n2 zeros(size(Xvalues{u},1),n3) -Xvalues{u} Xvalues{u} zeros(size(Xvalues{u}))];
    names{5*numcomp+u}      = ['Sa-Sp ' Xnames{u}];
end

% synesthetes proj - controls
for u=1:numcomp
    values{6*numcomp+u}     = [ones(size(Xvalues{u},1),n1)/n1 zeros(size(Xvalues{u},1),n2) -ones(size(Xvalues{u},1),n3)/n3 Xvalues{u} zeros(size(Xvalues{u})) -Xvalues{u}];
    names{6*numcomp+u}      = ['Sp-C ' Xnames{u}];
end

% synesthetes asso - controls
for u=1:numcomp
    values{7*numcomp+u}     = [zeros(size(Xvalues{u},1),n1) ones(size(Xvalues{u},1),n2)/n2 -ones(size(Xvalues{u},1),n3)/n3 zeros(size(Xvalues{u})) Xvalues{u} -Xvalues{u}];    
    names{7*numcomp+u}      = ['Sa-C ' Xnames{u}];
end

% synesthetes - controls
for u=1:numcomp
    values{8*numcomp+u}     = [ones(size(Xvalues{u},1),n1)/n1 ones(size(Xvalues{u},1),n2)/n2 2*(-ones(size(Xvalues{u},1),n3)/n3) Xvalues{u} Xvalues{u} -(2*Xvalues{u})];    
    names{8*numcomp+u}      = ['S-C ' Xnames{u}];
end

for n = 1 : numcomp-1
    types{n} = 'T';
    types{numcomp+n} = 'T';
    types{2*numcomp+n} = 'T';
    types{3*numcomp+n} = 'T';
    types{4*numcomp+n} = 'T';
    types{5*numcomp+n} = 'T';
    types{6*numcomp+n} = 'T';
    types{7*numcomp+n} = 'T';
    types{8*numcomp+n} = 'T';
end
for n = numcomp
    types{n}='F';
    types{numcomp+n} = 'F';
    types{2*numcomp+n} = 'F';
    types{3*numcomp+n} = 'F';
    types{4*numcomp+n} = 'F';
    types{5*numcomp+n} = 'F';
    types{6*numcomp+n} = 'F';
    types{7*numcomp+n} = 'F';
    types{8*numcomp+n} = 'F';
end

%%

cd(res_dir)
load SPM;
SPM = rmfield (SPM,'xCon');
% SPM.xCon = struct('name', {}, 'STAT', {}, 'c', {}, 'XO', {}, 'iXO', {}, 'X1o', {}, 'eidf', {}, 'Vcon', {}, 'Vspm', {});
for n=1:max(size(names))
    contrast = spm_FcUtil('Set',names{n}, types{n}, 'c', values{n}', SPM.xX.xKXs); % if error, check "cont" variable line 63, if good number of contrasts (including the resting when appropriate)
    SPM.xCon(n,1) = contrast;
end

save SPM SPM;

spm_contrasts(SPM);