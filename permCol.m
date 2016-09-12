function [perCol]= permCol(P)
%this function is to transform the permutation matrix into permutation
%vector. 
[mm,nn]=size(P);
perCol=[];
for i=1:mm
    for j=1:nn
        if P(i,j) == 1
            perCol=[perCol,j];
        end
    end
end
end

            