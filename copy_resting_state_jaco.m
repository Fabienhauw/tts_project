%% Copy paste anat files for Romain

from_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
to_dir   = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/rs_for_jaco';

cd(from_dir)
S = dir;
mask = ismember({S.name}, {'.', '..','meinfo.mat','list_subj'});
S(mask) = [];

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];

S = [S_droit ; S_con_app];

for k = 1 : length(S)
    directory = fullfile(from_dir, S(k).name);
    cd(directory)
    rs_file = 'RS/swf/s5wts_OC.nii'; %if isempty(anat), anat = dir('*T1*p2'); end; if length(anat)>1, anat = anat(end);, end; anat = anat(1).name;
    
    if ~exist (fullfile(to_dir, S(k).name, 'rs'))
        mkdir(fullfile(to_dir, S(k).name, 'rs'))
    end
    
    cd(directory)
    copyfile (rs_file, fullfile(to_dir,S(k).name, 'rs'))
end