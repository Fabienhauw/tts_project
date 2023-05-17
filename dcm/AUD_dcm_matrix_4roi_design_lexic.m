% to create all matrices possible for your DCM model;
nbreg=4;
lstgm=1; smg=2; sts=3; vwfa=4;

% basic matrices
A=eye(nbreg);
A(smg,sts) = 1;
A(sts,smg) = 1;
A(vwfa,lstgm) = 0;
A(vwfa,lstgm) = 0;
A(lstgm,vwfa) = 0; % no communication between the input and output rois;
A(lstgm,smg)  = 0;
A(lstgm,sts)  = 0; % the input has no feedback;
A(smg,vwfa)  = 0;
A(sts,vwfa)  = 0; % the output does not give feedback;

B=zeros(nbreg);
fix_conn_matrix = {
    smg,sts;
    sts,smg;
    vwfa,lstgm;
    lstgm,vwfa;
    lstgm,smg;
    lstgm,sts;
    smg,vwfa;
    sts,vwfa;
    lstgm,lstgm;
    smg,smg;
    sts,sts;
    vwfa,vwfa;
    };

nb_fix_conn1=size(fix_conn_matrix,1);
% for fixco=1:nb_fix_conn
%     A{1,16}(fix_conn{fixco,1},fix_conn{fixco,2});
% end


%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % you are at the case (row,col) of your matrix
        counter=0;
        for fixconn=1:nb_fix_conn1 % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix{fixconn,2}])
                counter=counter+1; %if the counter reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if counter==nb_fix_conn1
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
    if (a{1,aa}(smg,lstgm)==0 & a{1,aa}(sts,lstgm)==0)  | (a{1,aa}(vwfa,smg)==0 & a{1,aa}(vwfa,sts)==0)
        a{1,aa}=[];
    end
end
a=a(~cellfun('isempty',a));


%% matrices for modulation connectivity
lstgm=1; smg=2; sts=3; vwfa=4;
A=zeros(nbreg);

B1=zeros(nbreg);
B2=zeros(nbreg);

% version as in FLC, only connexions with STG are authorized
fix_conn_matrix1 = {
    lstgm,lstgm;
    smg,smg;
    sts,sts;
    vwfa,vwfa;
    vwfa,lstgm;
    vwfa, sts;
    vwfa, smg;
    lstgm,vwfa;
    lstgm,smg;
    lstgm,sts;
    sts, smg;
    smg, sts;
    smg,vwfa;
    sts,vwfa;
    };

% version with only connexions involving smg are authorized
% fix_conn_matrix1 = {
%     lstgm,lstgm;
%     smg,smg;
%     sts,sts;
%     vwfa,vwfa;
%     vwfa,lstgm;
%     vwfa, sts;
%     lstgm,vwfa;
%     lstgm,smg;
%     lstgm,sts;
%     smg,vwfa;
%     sts,vwfa;
%     sts,lstgm;
%     };

% first is for phonologic modulation: only connexions involving smg.

% version as in FLC, all connexions not involving stg are authorized
fix_conn_matrix2 = {
    lstgm,lstgm;
    smg,smg;
    sts,sts;
    vwfa,vwfa;% fin diagonale
    vwfa,lstgm;
    lstgm,vwfa;
    lstgm,smg;
    lstgm,sts;
    sts,lstgm;
    smg,vwfa;
    sts,vwfa;
    smg,lstgm;
    };


% version with only connexions involving sts are authorized
% fix_conn_matrix2 = {
%     lstgm,lstgm;
%     smg,smg;
%     sts,sts;
%     vwfa,vwfa;
%     vwfa,lstgm;
%     vwfa, smg;
%     lstgm,vwfa;
%     lstgm,smg;
%     lstgm,sts;
%     smg,vwfa;
%     sts,vwfa;
%     smg,lstgm;
%     };
% second is for lexical modulation: only connexions involving sts.

nb_fix_conn1=size(fix_conn_matrix1,1);
nb_fix_conn2=size(fix_conn_matrix2,1);

%% this part is to determine what are the cases that are authorized to vary their values;
for row=1:nbreg
    for col=1:nbreg % at this step, you are at the case (row,col) of your matrix
        counter=0;
        for fixconn=1:nb_fix_conn1 % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix1{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix1{fixconn,2}])
                counter=counter+1; %if the counter reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if counter==nb_fix_conn1
                    B1(row,col)=1;
                end
            else
                break
            end
        end
    end
end

for row=1:nbreg
    for col=1:nbreg % at this step, you are at the case (row,col) of your matrix
        counter=0;
        for fixconn=1:nb_fix_conn2 % then you check if this location is not a fixed connexion:
            if ~isequal([row], [fix_conn_matrix2{fixconn,1}]) | ~isequal ([col], [fix_conn_matrix2{fixconn,2}])
                counter=counter+1; %if the counter reaches the total number of fixed connexions, that means this location is different from all the fixed connexions locations;
                if counter==nb_fix_conn2
                    B2(row,col)=1;
                end
            else
                break
            end
        end
    end
end

%%
n = size(B1,1); m = size(B1,2);
for struct_model=1:length(a)
    constraints_model = (a{1,struct_model} == 1) & (B1 == 1); constraints_model = reshape(constraints_model', 1, n*m);
    s = 2^sum(constraints_model);
    result1 = repmat(constraints_model,s,1);
    result1(:,logical(constraints_model)) = dec2bin(0:s-1)-'0';
    
    constraints_model = (a{1,struct_model} == 1) & (B2 == 1); constraints_model = reshape(constraints_model', 1, n*m);
    s = 2^sum(constraints_model);
    result2 = repmat(constraints_model,s,1);
    result2(:,logical(constraints_model)) = dec2bin(0:s-1)-'0';
    
    cout_model = 0;
    
    for modul_model1 = 1 : size(result1,1)
        for modul_model2 = 1 : size(result2,1)
            cout_model = cout_model + 1 ;
            a_bis{cout_model,struct_model}(:,:,1) = a{1,struct_model};
            a_bis{cout_model,struct_model}(:,:,2) = reshape(result1(modul_model1,:), n, m)';
            a_bis{cout_model,struct_model}(:,:,3) = reshape(result2(modul_model2,:), n, m)';
        end
    end
end
clear a

a = a_bis;

nb_models = length(find(~cellfun(@isempty,a))); % this computes the total number of models;