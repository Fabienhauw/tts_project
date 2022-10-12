%% Mask
clear
clc

tts_group.path_to_subject;
tts_group.subj_choice;

for k = a:b
    cd (fullfile(D,S(k).name,'anat'))
    warpanat = dir('wmv*');
    warpanat = warpanat.name;
    fprintf('%s%s \n', 'MASK JOB ', S(k).name);
    unix(sprintf('%s%s','fsl5.0-bet ', warpanat, ' brain_extraction -R -f 0.3 -m'));
    unix(sprintf('%s%s', 'gunzip -f brain_extraction_mask.nii.gz'));
end