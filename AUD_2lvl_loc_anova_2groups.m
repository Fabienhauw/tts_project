clear;
clc;

% for any question, look at the function spm_run_factorial_design

spm('defaults','FMRI')
global defaults
global UFp; UFp = 0.001;
% categ='controls';
D = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images';
S = dir(D);
mask = ismember({S.name}, {'.', '..'});
S(mask) = [];

do_cov = 1;

Syn = S(~cellfun(@isempty,(regexp({S.name},'Sujet')))); 
Con = S(~cellfun(@isempty,(regexp({S.name},'Control'))));
S = [Syn;Con];

% left handed syn: Sujet05|Sujet07|Sujet11|Sujet14|Sujet16
% matched controls: Control02|Control04|Control05|Control07|Control17 

gaucher_appar = {'Control02|Control04|Control07|Control17|Control22|Control23|Control24|Control25|Control26|Sujet'};
% gaucher_appar = {'Control02|Control04|Control05|Control07|Control17|Sujet'};
mask_gauch_con = ~cellfun(@isempty,(regexp({S.name},gaucher_appar)));
S_con_app = S;
S_con_app(mask_gauch_con) = [];
con_group = repmat({'control'}, 1, length(S_con_app));

mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Sujet05|Sujet07|Sujet11|Sujet14|Sujet16|Control')));
% mask_gauch =  ~cellfun(@isempty,(regexp({S.name},'Control')));
S_droit = S;
S_droit(mask_gauch) = [];
syn_group = repmat({'syn'}, 1, length(S_droit));

categ = [syn_group, con_group];
nsub(1) = length(S_droit); % synesthetes
nsub(2) = length(S_con_app); % controls
ngroups = 2;

S_effect = [S_droit ; S_con_app];

if do_cov
    res_dir_base = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ANOVA_s5_without_resting_with_pcs');
else
    res_dir_base = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ANOVA_s5_without_resting_%s_to_%s_s8', S_effect(1).name, S_effect(end).name);
end

% res_dir_base = sprintf('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/ANOVA_s5_without_resting_test');


vector_age = [
    25.1013699; 70.8219178; 23.6821918; 24.3342466; 21.109589; 31.7260274; 18.8438356; ...
    36.3424658; 49.7753425; 27.2767123; 26.3616438; 40.8876712; 22.8246575; 43.1013699; 44.309589; ...
    51.2246575; 40.7561644; 18.5835616; 39.939726; 58.9534247; 43.1753425; 39.9260274; ... % end of synesthetes
    35.8219178; 23.2054795; 21.9726027; 31.7945205; 30.7589041; 70.3808219; 26.3589041; ... %start of controls
    22.660274; 42.0027397; 45.5589041; 19.6027397; 55.9041096; 19.4958904; 38.3643836; ...
    49.9178082; 46.0739726; 51.8630137; 25.1945205; 23.0547945; 41.4383562; 41.5972603; ...
    29.865753; 26.57534247; 24.67945205; 29.72328767; 58.97534247; ...
    ];

vector_hand = [
    0; 0; 0; 0; 1; 0; 1; 0; 0; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0; 0; 0;... % end of synesthetes
    zeros(21,1); 1; 1; 1; 1; 1; ... % end of controls
    ]; %0 = right, 1 = left;

vector_pcs = [
    0;0;1;1;1;1;1;1;0;0;1;1;0;1;1;0;1;1;0;0;0;1;1;0;0;0;0;0;0;1;0;0;0;1; ...
    ];

for j = 1 : size(S,1)
    if ~isempty(find(~cellfun(@isempty,(regexp({S_effect.name},S(j).name)))))
        mask_cov(j,1) = 1;
    else
        mask_cov(j,1) = 0;
    end
end

% vector_cov1 = vector_age(mask_cov==1); name_cov1 = 'age';
vector_cov1 = vector_pcs; name_cov1 = 'pcs';
vector_cov2 = vector_hand(mask_cov==1);
subname = {S_effect.name};

totsub=length(subname);

%-----------------------------------------------------------------
cd(fullfile(D,S(1).name, 'Aud/loc/stats_s5_without_resting'));
cont = [1 : 5];
% words, pseudowords, numbers, normal_speech, scramble_speech, odds, motor
%-----------------------------------------------------------------

ncon = length(cont);
nscan = ncon*totsub;
if ~isdir(res_dir_base)
    mkdir(res_dir_base)
end

vector_cov1 = repmat(vector_cov1, 1, ncon);
vector_cov1 = reshape(vector_cov1, size(vector_cov1,1)*size(vector_cov1, 2), 1);

cd(res_dir_base)

% get image files names
P={};
for con=1:ncon
    for s=1:totsub
        sub=subname{s};
        P{(con-1)*totsub+s} =	sprintf(['/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/final_images/%s/Aud/loc/stats_s5_without_resting/s8con_%04d.nii'],sub,cont(con));
    end
end

j=0;
for i=1:length(P)
    if ~exist(P{i})
        j=j+1;
        filestosmooth{j}=strrep(P{i},'s8con','con');
    end
end

if j>0
    for u=1:j
        spm_smooth(filestosmooth{u},strrep(filestosmooth{u},'con_','s8con_'),[8 8 8],0); 
    end
end

%% -Assemble SPM structure
%=======================================================================

SPM.nscan = nscan;
SPM.xY.P = P;

for i=1:SPM.nscan
    SPM.xY.VY(i) = spm_vol(SPM.xY.P{i});
end

cname = cell(ncon,1);
k=0;
for c=1:ncon
    for group=1:ngroups
        k=k+1;
        cname{k} = sprintf('con %d group %d',c,group);
    end
end

if do_cov == 1
    k=k+1;
    cname{k} = 'age';
end

for group=1:ngroups
    for c=1:nsub(group)
        k=k+1;
        cname{k} = sprintf('sub %d group %d',c,group);
    end
end

% for c=1:ncon
%     for group=1:ngroups
%         for tmp_sub=1:nsub(group)
%             k=k+1;
%             cname{k} = sprintf('con %d sub %d group %d', c, tmp_sub, group);
%         end
%     end
% end

%%%% define the subject part of the matrix
oc = ones(ncon,1);
subjectpart=kron(oc,eye(totsub));
subjectpart=subjectpart - repmat(mean(subjectpart),ncon*totsub,1);

%%% define the contrast part of the matrix
os = zeros(totsub,ngroups);
g=0;
for group=1:ngroups
    os((g+1):(g+nsub(group)),group)=1;
    g=g+nsub(group);
end

% SPM.xX = struct(...
%     'X',[kron(eye(ncon),os) subjectpart] ,... % change this
%     'iH',1:ncon,'iC',zeros(1,0),...
%     'iB',ncon+[1:totsub],...
%     'iG',zeros(1,0),...
%     'name',{cname},'I',[ones(nscan,1) kron([1:ncon]',os) kron(oc,[1:totsub]') ones(nscan,1)],...
%     'sF',{{'repl'  'cond'  'subj'  ''}});


rname_cov1 = name_cov1;
rvector_cov1 = vector_cov1;

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
%--------------------------------------------------------------------------
sCC = {'around overall mean';...                            %-1
    'around sF1 means';...                                  %-2
    'around sF2 means';...                                  %-3
    'around sF3 means';...                                  %-4
    'around sF4 means';...                                  %-5
    'around sF2 (within sF4) means';...                     %-6
    'around sF3 (within sF4) means';...                     %-7
    '<no centering>';...                                    %-8
    'around user specified value';...                       %-9
    '(as implied by AnCova)';...                            %-10
    'GM';...                                                %-11
    '(redundant: not doing AnCova)'}';                      %-12

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
%--------------------------------------------------------------------------
CFIforms = {'[]',   'C',    '{}';...                        %-1
    'I(:,1)',       'FxC',  '{sF{1}}';...                   %-2
    'I(:,2)',       'FxC',  '{sF{2}}';...                   %-3
    'I(:,3)',       'FxC',  '{sF{3}}';...                   %-4
    'I(:,4)',       'FxC',  '{sF{4}}';...                   %-5
    'I(:,[4,2])',   'FxC',  '{sF{4},sF{2}}';...             %-6
    'I(:,[4,3])',   'FxC',  '{sF{4},sF{3}}' };              %-7

CCforms = {'ones(nscan,1)',CFIforms{2:end,1},''}';
vector_cov1 = rvector_cov1 - spm_meanby(rvector_cov1,eval(CCforms{1}));

name_cov1     = {name_cov1};
iCC = 1;
iCFI = 1;

str = {sprintf('%s',rname_cov1)};
if size(rvector_cov1,2)>1, str = {sprintf('%s (block of %d covariates)',...
        str{:},size(rvector_cov1,2))}; end
if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

typ = 1;
if do_cov == 1
    SPM.xC = struct(...
        'rc', rvector_cov1,         'rcname', rname_cov1, ...
        'c', vector_cov1,           'cname', {name_cov1}, ...
        'iCC', 1,         'iCFI', 1, ...
        'type', typ, ...
        'cols', ncon*size(os,2)+1,... %[1:size(vector_cov1,2)] + size([H,C],2) + size([B,G],2)*min(typ-1,1), ...
        'descrip', {{'age'; 'used centered around overall mean'}});
else
    SPM.xC = [];
end

if do_cov == 1
    SPM.xX = struct(...
        'X',[kron(eye(ncon),os) vector_cov1 subjectpart] ,... % change this
        'iH',1:ncon, ...
        'iC', ncon+1, ...
        'iB',ncon+1+[1:totsub],...
        'iG',zeros(1,0),...
        'name',{cname},'I',[ones(nscan,1) kron([1:ncon]',os) kron(oc,[1:totsub]') ones(nscan,1)],...
        'sF',{{'repl'  'cond'  'subj'  ''}});
    
else
    SPM.xX = struct(...
        'X',[kron(eye(ncon),os) subjectpart] ,... % change this
        'iH',1:ncon,'iC',zeros(1,0),...
        'iB',ncon+[1:totsub],...
        'iG',zeros(1,0),...
        'name',{cname},'I',[ones(nscan,1) kron([1:ncon]',os) kron(oc,[1:totsub]') ones(nscan,1)],...
        'sF',{{'repl'  'cond'  'subj'  ''}});
    
end

SPM.xGX = struct(...
    'iGXcalc',1,    'sGXcalc','omit',                               'rg',[],...
    'iGMsca',9,     'sGMsca','<no grand Mean scaling>',...
    'GM',0,         'gSF',ones(nscan,1),...
    'iGC',  12,     'sGC',  '(redundant: not doing AnCova)',        'gc',[],...
    'iGloNorm',9,   'sGloNorm','<no global normalisation>');
donotuseVI = 1;
if (donotuseVI)
    SPM.xVi = struct(...
        'iid',1, 'V',speye(ncon*totsub) );
else
    x=zeros(ncon);
    s=eye(totsub);
    nv=0;
    % first do the leading diagonal elements
    for group=1:ngroups
        for i=1:ncon
            for j=i:ncon
                nv=nv+1;
                %v=x;
                %v(j,i)=1; v(i,j)=1;
                vitemp=zeros(ncon*totsub,ncon*totsub);
                deci = (i-1)*totsub;
                decj = (j-1)*totsub;
                for g2=1:group-1
                    deci =deci+nsub(g2);
                    decj =decj+nsub(g2);
                end
                for s=1:nsub(group)
                    vitemp(  deci+s,  decj+s)=1;
                    vitemp(  decj+s,  deci+s)=1;
                end
                vi{nv}=vitemp;
            end
        end
    end
    SPM.xVi = struct(...
        'iid',0, 'I',SPM.xX.I, 'sF','SPM.xX.sF',...
        'var',[0 1 0 0], 'dep',[0 1 0 0], 'Vi',{vi} );
end

%- With IMPLICIT masking
%=======================================================================
%     Mdes 	= struct(	'Analysis_threshold',	{'None (-Inf)'},...
%         'Implicit_masking',	{'Yes: NaNs treated as missing'},...
%         'Explicit_masking',	{'No'});
%     SPM.xM	= struct(	'T',-Inf,'TH',ones(nsub,1)*-Inf,...
%         'I',1,'VM',[],'xs',Mdes);

%- With EXPLICIT masking
%=======================================================================
Mdes 	= struct(	'Analysis_threshold',	{'None (-Inf)'},...
    'Implicit_masking',	{'No'},...
    'Explicit_masking',	{'Yes'});
exp_mask = fullfile('/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/masks',S_effect(1).name);
exp_mask = sprintf('%s_%s_aud_loc_mask_thr_s5.nii',exp_mask, S_effect(end).name);
SPM.xM	= struct(	'T',-Inf,'TH',ones(nscan,1)*-Inf,...
    'I',1,'VM',spm_vol([exp_mask]),'xs',Mdes);

Pdes    = {{sprintf('%d condition, +1  covariate, +0 block, +0 nuisance',ncon); sprintf('%d total, having %d degrees of freedom',ncon,ncon); sprintf('leaving %d degrees of freedom from %d images',nscan-ncon,nscan)}};

SPM.xsDes = struct(...
    'Design',               {'1-way ANOVA (within-between-subjects)'},...
    'Global_calculation',   {'omit'},...
    'Grand_mean_scaling',   {'<no grand Mean scaling>'},...
    'Global_normalisation', {'<no global normalisation>'},...
    'Parameters',           Pdes);

SPM.SPMid       = 'SPM12: spm_spm (v7738)';

% spm_run_factorial_design(SPM)
% spm_design_within_subject(fblock,cov)

save SPM SPM


%%
% Estimate parameters
%===========================================================================
SPM = spm_spm(SPM);

%%
%%%%%% Second, define the contrasts and estimate them

words           = [1 0 0 0 0];
pseudowords     = [0 1 0 0 0];
numbers         = [0 0 1 0 0];
normal_speech   = [0 0 0 1 0];
scramble_speech = [0 0 0 0 1];
% odds            = [0 0 0 0 0 1 0];
% motor           = [0 0 0 0 0 0 1];
EOI             = [eye(5)];

% words           = [1 zeros(1,ncon-1)];
% pseudowords     = zeros(1,length(words)); pseudowords(2) = 1;
% numbers         = zeros(1,length(words)); numbers(3) = 1;
% normal_speech   = zeros(1,length(words)); normal_speech(4) = 1;
% scramble_speech = zeros(1,length(words)); scramble_speech(5) = 1;
% odds            = zeros(1,length(words)); odds(6) = 1;
% motor           = zeros(1,length(words)); motor(7) = 1;

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

% % expand contrasts for a model with derivatives:
% [a,b]=size(parameters.contrast.values);
% for u=1:length(parameters.contrast.values)
% parameters.contrast.values{u}=reshape([parameters.contrast.values{u};zeros(1,b)],1,2*b);
% end

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
    values{u}               = reshape( [Xvalues{u};zeros(size(Xvalues{u}))],     size(Xvalues{u},1) , 2*size(Xvalues{u},2));
%     values{u}               = reshape( [repmat(Xvalues{u}, nsub(1),1); repmat(zeros(size(Xvalues{u})), nsub(2),1)],     size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{u}                = ['S ' Xnames{u}];
end
% controls
for u=1:numcomp
    values{numcomp+u}       = reshape([zeros(size(Xvalues{u})); Xvalues{u}],size(Xvalues{u},1),2*size(Xvalues{u},2));
%     values{numcomp+u}       = reshape([repmat(zeros(size(Xvalues{u})), nsub(1),1); repmat(Xvalues{u}, nsub(2),1)],     size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{numcomp+u}        = ['C ' Xnames{u}];
end
% (synesthetes+controls)
for u=1:numcomp
    values{2*numcomp+u}     = reshape([Xvalues{u};Xvalues{u}],size(Xvalues{u},1),2*size(Xvalues{u},2));
%     values{2*numcomp+u}     = reshape([repmat(Xvalues{u}, nsub(1),1); repmat(Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{2*numcomp+u}      = ['S+C ' Xnames{u}];
end
% synesthetes-controls
for u=1:numcomp
    values{3*numcomp+u}     = reshape([Xvalues{u};-Xvalues{u}],size(Xvalues{u},1),2*size(Xvalues{u},2));
%     values{3*numcomp+u}     = reshape([repmat(Xvalues{u}, nsub(1),1); repmat(-Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{3*numcomp+u}      = ['S-C ' Xnames{u}];
end
% controls-synesthetes
for u=1:numcomp
    values{4*numcomp+u}     = reshape([-Xvalues{u};Xvalues{u}],size(Xvalues{u},1),2*size(Xvalues{u},2));    
%     values{4*numcomp+u}     = reshape([repmat(-Xvalues{u}, nsub(1),1); repmat(Xvalues{u}, nsub(2),1)], size(Xvalues{u},1) , totsub*size(Xvalues{u},2));
    names{4*numcomp+u}      = ['C-S ' Xnames{u}];
end

if do_cov == 1 
    for n=1:size(values,2)
        values{n}=[values{n} zeros(size(values{n},1), size(name_cov1,1)) zeros(size(values{n},1),totsub)];
    end
    % one contrast for each cov
    n = n + 1;
    values{n}=[zeros(size(Xvalues{1})) zeros(size(Xvalues{1})) 1 zeros(1,totsub)];
    names{n} = rname_cov1;
    
    n = n + 1;
    values{n}=[zeros(size(Xvalues{1})) zeros(size(Xvalues{1})) -1 zeros(1,totsub)];
    names{n} = ['-' rname_cov1];
    
else
    for n=1:size(values,2)
        values{n}=[values{n} zeros(size(values{n},1),totsub)];
    end
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

if do_cov == 1
    types{5*numcomp+1} = 'T';
    types{5*numcomp+2} = 'T';
end

%%

load SPM;
SPM = rmfield (SPM,'xCon');
% SPM.xCon = struct('name', {}, 'STAT', {}, 'c', {}, 'XO', {}, 'iXO', {}, 'X1o', {}, 'eidf', {}, 'Vcon', {}, 'Vspm', {});
for n=1:max(size(names))
    contrast = spm_FcUtil('Set',names{n}, types{n}, 'c', values{n}', SPM.xX.xKXs); % if error, check "cont" variable line 63, if good number of contrasts (including the resting when appropriate)
    SPM.xCon(n,1) = contrast;
end

save SPM SPM;

spm_contrasts(SPM);