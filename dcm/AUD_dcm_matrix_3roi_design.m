% to create all matrices possible for your DCM model;
nbreg=3;
lpstg=1; smg=2; vwfa=3; 

% basic matrices
A=eye(nbreg);
A(smg,vwfa) = 1;
A(vwfa,smg) = 1;
A(lpstg,smg) = 1;
A(smg,lpstg) = 1;
A(vwfa,lpstg) = 0;
A(lpstg,vwfa) = 0; % stg never receives feedback;

B=zeros(nbreg);
fix_conn_matrix = {
    smg,vwfa;
    vwfa,smg;
    lpstg,smg;
    smg,lpstg;
    vwfa,lpstg;
    lpstg,vwfa;
    lpstg,lpstg;
    smg,smg;
    vwfa,vwfa;
    };

nb_fix_conn=size(fix_conn_matrix,1);
% for fixco=1:nb_fix_conn
%     A{1,16}(fix_conn{fixco,1},fix_conn{fixco,2});
% end


%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % you are at the case (row,col) of your matrix
        counter=0;
        for fixconn=1:nb_fix_conn % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix{fixconn,2}])
                counter=counter+1; %if the counter reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if counter==nb_fix_conn
                    B(row,col)=1;
                end
            else
                break
            end
        end
    end
end

%%
count=1;
aa=A;
a{count}=aa;

nb_struc_mod=length(a);

% now you have to delete models with no connexions from the stgm or no
% connexions to the vwfa:
for aa=1:nb_struc_mod
    if a{1,aa}(smg,lpstg)==0 & a{1,aa}(vwfa,lpstg)==0
        a{1,aa}=[];
    end
end
a=a(~cellfun('isempty',a));

%% matrices for modulation connectivity
lpstg=1; smg=2; vwfa=3;
A=zeros(nbreg);

B=zeros(nbreg);
fix_conn_matrix = {
    lpstg,lpstg;
    smg,smg; % mtg and smg are always connected
    vwfa,vwfa;
    vwfa,lpstg;
    lpstg,vwfa;
    };

nb_fix_conn=size(fix_conn_matrix,1);

%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % at this step, you are at the case (row,col) of your matrix
        counter=0;
        for fixconn=1:nb_fix_conn % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix{fixconn,2}])
                counter=counter+1; %if the counter reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if counter==nb_fix_conn
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
                    con_count = con_count -1;
                end
            end
            con_count = con_count -1;
        end
    end
end

nb_models = length(find(~cellfun(@isempty,a))); % this computes the total number of models;