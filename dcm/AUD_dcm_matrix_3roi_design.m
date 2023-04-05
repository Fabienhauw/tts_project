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
n = size(B,1); m = size(B,2);

constraints_model = reshape(B', 1, n*m);
s = 2^sum(constraints_model);
result = repmat(constraints_model,s,1);
result(:,logical(constraints_model)) = dec2bin(0:s-1)-'0';
for struct_model = 1 : size(result,1)
    a{1,struct_model} = reshape(result(struct_model,:), n, m)';
end

for tmp_a = 1 : length(a)
    a{tmp_a} = a{tmp_a} + A;
end

for aa=1:length(a)
    if (a{1,aa}(smg,lpstg)==0 & a{1,aa}(mfg,lpstg)==0)  | (a{1,aa}(vwfa,smg)==0 & a{1,aa}(vwfa,mfg)==0)
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
n = size(B,1); m = size(B,2);
for struct_model=1:length(a)
    constraints_model = (a{1,struct_model} == 1) & (B == 1); constraints_model = reshape(constraints_model', 1, n*m);
    s = 2^sum(constraints_model);
    result = repmat(constraints_model,s,1);
    result(:,logical(constraints_model)) = dec2bin(0:s-1)-'0';
    for modul_model = 1 : size(result,1)
        a{modul_model + 1,struct_model} = reshape(result(modul_model,:), n, m)';
    end
end

nb_models = length(find(~cellfun(@isempty,a))); % this computes the total number of models;