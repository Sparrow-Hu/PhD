function [ Conflicting,Redundant ] = luRC( table,tol )
%LURC use LU factorization to find conflicting and redundant equations
%   LU factorization is an variant of Gauss Jordan Elimination 
Conflicting={};
Redundant={};
matrix=table2array(table(:,1:end-1));
r=rank(matrix);
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

end

