% % Example of function getsubjects.data.m :

%%%% give the number of subject groups
ngroups=2;
%%%% name the groups of subjects 
groupn{1} = 'Synesthetes';
groupn{2} = 'Controls';
%%%% give the number of subjects in each group
nsub(1) = 17;
nsub(2) = 17;
totsub = sum(nsub);  %% total number of subjects
res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/FIR_analyses_rh';
%% first in simplified form:
D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
mask = ismember({S.name}, {'.', '..'});
S(mask) = [];
gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
mask_gauch =  ~cellfun(@isempty,(regexp({S.name},{'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control'})));
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];

S_droit = S;
S_droit(mask_gauch) = [];

S = [S_droit ; S_con_app];
subjects = {S.name};
%% then with the full path
modeldir = 'Aud/loc/stats_s5_without_resting';
for k=1:totsub
    subjectsdir{k} = fullfile(D, S(k).name, modeldir);
end