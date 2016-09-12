function [Conflicting,Redundant] = findRC( table,tol )
%FINDRC Find the number of overconstraints as well as output the distance
%between b and B
%   SVD---find the distribution of singular values
%   QR---find the ordering of equations in accordance with the distribution
%   of singular values
%   LU--Find the linear relationship between equations;this LU
%   decomposition is based on the ordering given by QR factorization
Conflicting={};
Redundant={};
matrix1=table2array(table(:,1:end-1));
mat=matrix1';
%% SVD to find the distribution of singular values: the value for the gap is : 1e5
S=svd(matrix1);
r=rank(matrix1);
semilogy(diag(S),'r.');
l=length(S);
trucat=[];
% for i=1:l-1
%     if S(i)/S(i+1)>=1e3
%         trucat = [trucat,i];
%     end   
% end
%perhaps a better way is to cluster these singular values
%% equation order after QR factorization
[Q,R,p]=qr(mat,0); 
if isempty(trucat)
    table=table(p,:);
    matrix=table2array(table(:,1:end-1));
    [L,U,P]=lu(matrix);
    P=permCol(P);
    table=table(P,:);
    [m,n]=size(L);
    b=table{:,end};
    for i=1:n % starting from the column
        if i+1<=m
            for j=i+1:m %iterate along row direction
                b(j)=b(j)-L(j,i)*b(i);
            end
        end
    end
    table{:,end}=b;
    if r+1<=m
        temp1=1;
        temp2=1;
        for i=r+1:m
            if abs(b(i)) >=tol
                Conflicting(temp1) = table.Properties.RowNames(i);
                 temp1=temp1+1;
            else
                Redundant(temp2) = table.Properties.RowNames(i);
                temp2=temp2+1;
            end
        end
    end
    
else
table0=table(p(1:trucat(1)),:);% Basis constraints
matrix0=rref(table2array(table0));
table1=table(p(trucat(1)+1:trucat(2)),:);% Basis + ill constraints1
table2=table([p(1:trucat(1)),p(trucat(2)+1:l)],:); % Basis + ill constraints2
table3=table([p(1:trucat(1)),p(l+1:end)],:); % Basis + overconstraints
%% LU decomposition to find linear relationship between equations in table1
matrix1 = table2array(table1);
[L,U,P]=lu([matrix0;matrix1]);
P=permCol(P);%transform the permutation matrix to permutation vector
table1=table1(P,:);%reorder the rows based on permutation matrix given by LU factorization
matrix1=table2array(table1);
S2=svd(matrix1);
%% LU decomposition to find linear relationship between equations in table2
matrix2 = table2array(table2);
[L,U,P]=lu(matrix2);
P=permCol(P);%transform the permutation matrix to permutation vector
table2=table2(P,:);%reorder the rows based on permutation matrix given by LU factorization
matrix2=table2array(table2);
S2=svd(matrix2);
%% LU decomposition to find linear relationship between equations in table3
matrix3 = table2array(table3);
[L,U,P]=lu(matrix3);
P=permCol(P);%transform the permutation matrix to permutation vector
table3=table3(P,:);%reorder the rows based on permutation matrix given by LU factorization
matrix3=table2array(table3);
S2=svd(matrix3);
end