% to create all matrices possible for your DCM model;
% mtg towards vwfa and smg towards vwfa, mtg towards smg and smg towards
% mtg modulations are suppressed here: only modulations from stg towards
% mtg and smg are present.

nbreg=4;
stgm=1; smg=2; mtg=3; vwfa=4;

% basic matrices
A=eye(nbreg);
A(smg,mtg) = 1;
A(mtg,smg) = 1;
A(vwfa,stgm) = 0;
A(smg,vwfa) = 0;
A(mtg,vwfa) = 0; % vwfa nevers sends back connectivity;
A(stgm,vwfa) = 0;
A(stgm,smg) = 0;
A(stgm,mtg) = 0; % stgm never receives feedback;

B=zeros(nbreg);
fix_conn_matrix = {
    smg,mtg;
    mtg,smg; % mtg and smg are always connected
    vwfa,stgm;
    smg,vwfa;
    mtg,vwfa; % vwfa nevers sends back connectivity;
    stgm,vwfa;
    stgm,smg;
    stgm,mtg; % stgm never receives feedback;
    stgm,stgm;
    smg,smg;
    mtg,mtg;
    vwfa,vwfa;
    };

nb_fix_conn=size(fix_conn_matrix,1);
% for fixco=1:nb_fix_conn
%     A{1,16}(fix_conn{fixco,1},fix_conn{fixco,2});
% end


%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % at this step, you are at the case (row,col) of your matrix
        sum=0;
        for fixconn=1:nb_fix_conn % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix{fixconn,2}])
                sum=sum+1; %if the sum reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if sum==nb_fix_conn
                    B(row,col)=1;
                end
            else
                break
            end
        end
    end
end

%%
variabl_conn = find(B==1);
num_variabl_conn = length(variabl_conn);
count=1;
con_count=1;
for value=1:2
    aa=A;
    aa(variabl_conn(con_count))=value-1;
    if con_count<num_variabl_conn
        con_count = con_count + 1;
        for value=1:2
            aa(variabl_conn(con_count))=value-1;
            if con_count<num_variabl_conn
                con_count = con_count + 1;
                for value=1:2
                    aa(variabl_conn(con_count))=value-1;
                    if con_count<num_variabl_conn
                        con_count = con_count + 1;
                        for value=1:2
                            aa(variabl_conn(con_count))=value-1;
                            a{count}=aa;
                            count=count+1;
                        end
                        con_count = con_count -1;
                    end
                end
                con_count = con_count -1;
            end
        end
        con_count = con_count -1;
    end
end

nb_struc_mod=length(a);

% now you have to delete models with no connexions from the stgm or no
% connexions to the vwfa:
for aa=1:nb_struc_mod
    if a{1,aa}(smg,stgm)==0 & a{1,aa}(mtg,stgm)==0 | a{1,aa}(vwfa,smg)==0 & a{1,aa}(vwfa,mtg)==0
        a{1,aa}=[];
    end
end
a=a(~cellfun('isempty',a));

%% note : essayer avec boucle while ?

%% matrices for modulation connectivity
stgm=1; smg=2; mtg=3; vwfa=4;
A=zeros(nbreg);

B=zeros(nbreg);
fix_conn_matrix = {
    stgm,stgm;
    smg,smg; 
    mtg,mtg;
    vwfa,vwfa;
    stgm,smg; 
    stgm,mtg;
    stgm,vwfa;
    vwfa,stgm;
    mtg,vwfa; 
    smg,vwfa;
    vwfa,mtg;
    vwfa,smg;
    smg,mtg;
    mtg,smg;
    };

nb_fix_conn=size(fix_conn_matrix,1);

%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % at this step, you are at the case (row,col) of your matrix
        sum=0;
        for fixconn=1:nb_fix_conn % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix{fixconn,2}])
                sum=sum+1; %if the sum reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if sum==nb_fix_conn
                    B(row,col)=1;
                end
            else
                break
            end
        end
    end
end

%%
variabl_conn = find(B==1);
num_variabl_conn = length(variabl_conn);
for struct_model=1:length(a)
    count=1;
    con_count=1;
    aa=A;
    for value=1:2
        if a{1,struct_model}(variabl_conn(con_count))==0 & value==1
            aa(variabl_conn(con_count))=0;
        elseif a{1,struct_model}(variabl_conn(con_count))==0 & value==2 %before adding modulatory input on this
            % connexion, you have to be sure that this connexion exists in your structural model;
            % If the struc connexion does not exist, the only modulation
            % value for this connexion is 0, so no need for a
            % "second tour" and a lot a new identical models except for this connexion
            break
        else
            aa(variabl_conn(con_count))=value-1;
        end
        if con_count<num_variabl_conn
            con_count = con_count + 1;
            for value=1:2
                if a{1,struct_model}(variabl_conn(con_count))==0 & value==1
                    aa(variabl_conn(con_count))=0;
                elseif a{1,struct_model}(variabl_conn(con_count))==0 & value==2 %before adding modulatory input on this
                    % connexion, you have to be sure that this connexion exists in your structural model;
                    % If the struc connexion does not exist, the only modulation
                    % value for this connexion is 0, so no need for a
                    % "second tour" and a lot a new identical models except for this connexion
                    break
                else
                    aa(variabl_conn(con_count))=value-1;
                end
                a{1+count,struct_model}=aa;
                count=count+1;
            end
            con_count = con_count -1;
        end
    end
end

nb_models = length(find(~cellfun(@isempty,a))) - size(a,2); % this computes the total number of models;