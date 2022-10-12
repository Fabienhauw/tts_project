% https://github.com/CMRR-C2P/MB

clear
clc

addpath(genpath('/network/lustre/iss02/home/fabien.hauw/Documents/matvol'))
addpath(genpath('/network/lustre/iss01/home/fabien.hauw/Documents/MATLAB'));

main_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/NIFTI';
% subj_to_preproc = '((SYNESTHEX_Sujet)|(SYNESTHEX_Control))';
subj_to_preproc = 'Control18';
% subj_to_preproc = 'Sujet\d[4-9]|Sujet[1-2][0-9]';
e = exam(main_dir, subj_to_preproc); % all subjects with multi-echo

e.addSerie('PhysioLog$','physio')

e.getSerie('physio').addVolume('dcm$','dcm',1)

physio_file = e.getSerie('physio').getVolume('dcm').removeEmpty.getPath';

for f = 1 : length(physio_file)
    
    extractCMRRPhysio(physio_file{f})

end 
