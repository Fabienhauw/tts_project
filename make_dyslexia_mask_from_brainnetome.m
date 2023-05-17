%%
clear;clc;
matlabbatch{1}.spm.util.imcalc.input = {'/network/lustre/iss02/cohen/data/Fabien_official/atlas/brainnetome/BN_Atlas_246_1mm_reoriented.nii,1'};
matlabbatch{1}.spm.util.imcalc.output = 'dyslexia_network_brainnetome';
matlabbatch{1}.spm.util.imcalc.outdir = {'/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/masks'};
matlabbatch{1}.spm.util.imcalc.expression = '9*(i1==9) + 20*(i1==20) + 32*(i1==32) + 36*(i1==36) + 54*(i1==54) + 61*(i1==61) + 62*(i1==62) + 67*(i1==67) + 69*(i1==69) + 70*(i1==70) + 71*(i1==71) + 73*(i1==73)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', matlabbatch)
