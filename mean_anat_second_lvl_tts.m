tts_group.path_to_subject;
wd = pwd;

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
% mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},'Control02|Control04|Control07|Control17|Sujet')));
S_con_app = S;
S_con_app(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];

S_effect = [S_droit ; S_con_app];

clear matlabbatch
i=1;

experiment = {'vis_loc', 'aud_loc', 'vis_col'};
experiment_path = {'Vis/loc', 'Aud/loc', 'Vis/unfr_col'};
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Chose which comparaison you want %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = S_effect;

exp = 1;
exp_name = experiment{exp};
exp_path = experiment_path{exp};
%%
scans = {};
masks = {};
for k=1:length(S)
    cd(fullfile(D,S(k).name,'anat'))
    vol = dir ('wm*');
    vol = vol.name;
    vol = fullfile(D,S(k).name,'anat',vol);
    scans = [scans;vol];
    
    cd (fullfile(D,S(k).name, exp_path, '/stats_s5'))
    mask = dir('mask*');
    mask = mask.name;
    mask = fullfile(D,S(k).name,exp_path, '/stats_s5', mask);
    masks = [masks;mask];
end

outdir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/anat';
if ~isdir(outdir)
    mkdir(outdir)
end

matlabbatch{1}.spm.util.imcalc.input = scans;
matlabbatch{1}.spm.util.imcalc.output = sprintf('%s_%s_mean_anat',S(1).name, S(end).name);
matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
total = '';
for k = 1:length(S)
    map = sprintf('%s%d','i',k);
    if k == 1
        total = map;
    else
        total = [total '+' map];
    end
end
total = ['(' total ')'];
matlabbatch{1}.spm.util.imcalc.expression = [total '/' num2str(length(S))];
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

%%
matlabbatch{2}.spm.util.imcalc.input = masks;
matlabbatch{2}.spm.util.imcalc.output = sprintf('%s_%s_mask_s5',S(1).name, S(end).name);
if ~isdir('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks')
    mkdir('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks');
end
matlabbatch{2}.spm.util.imcalc.outdir = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks'};
total = '';
for k = 1:length(S)
    map = sprintf('%s%d','i',k);
    if k == 1
        total = map;
    else
        total = [total '+' map];
    end
end
total = ['(' total ')'];
matlabbatch{2}.spm.util.imcalc.expression = [total '/' num2str(length(S))];
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = -7;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;

%%
mask_input = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S(1).name);
mask_input = sprintf('%s_%s_mask_s5.nii',mask_input, S(end).name);
matlabbatch{3}.spm.util.imcalc.input = {mask_input};
matlabbatch{3}.spm.util.imcalc.output = sprintf('%s_%s_%s_mask_thr_s5',S(1).name, S(end).name, exp_name);
matlabbatch{3}.spm.util.imcalc.outdir = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks'};
matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.5';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = -7;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch)

cd(wd)
