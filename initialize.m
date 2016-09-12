function [ T ] = initialize( equation,variables,B )
%initialize: transforming the equations to adjacent table, where the last
%column is 'b'; column 1:end-1 is the adjacent matrix
%   first compute the adjacent matrix:admatrix
m=length(equation);
n=length(variables);
varCell={};
admatrix = zeros(m,n);
for i=1:m
    eq = equation(i);
    var_eq = symvar(eq);
    k =length(var_eq);
    for j=1:n
        for jj=1:k
            if strcmp(char(var_eq(jj)),char(variables(j)))
               admatrix(i,j) = 1;
            end
        end
    end
end
for i=1:n
    varCell=[varCell,char(variables(i))];
end
T =array2table(admatrix,'VariableNames',varCell); 
T.b=B';
%%second decompose the adjacent matrix using DM algorithm
end

