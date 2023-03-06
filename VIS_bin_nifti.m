% to binarize maps
clear
map_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/images/vis_anova';
cd(map_dir)
maps = dir('*.nii');
i=1;
for tmp_m = 1 : length(maps)
    cmap = maps(tmp_m).name;
    matlabbatch{i}.spm.util.imcalc.input = {fullfile(map_dir , cmap)};
    matlabbatch{i}.spm.util.imcalc.output = ['bin_' cmap];
    matlabbatch{i}.spm.util.imcalc.outdir = {map_dir};
    matlabbatch{i}.spm.util.imcalc.expression = 'i1>0';
    matlabbatch{i}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{i}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{i}.spm.util.imcalc.options.mask = 0;
    matlabbatch{i}.spm.util.imcalc.options.interp = 1;
    matlabbatch{i}.spm.util.imcalc.options.dtype = 4;
    i = i + 1;
end

spm_jobman('run', matlabbatch)