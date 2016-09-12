function [ adMatrix ] = adMatrix( matrix )
%ADMATRIX turn coefficient matrix into adjacent matrix
%   Detailed explanation goes here
[m,n]=size(matrix);
for i=1:m
    for j=1:n
        if matrix(i,j) ~= 0
            matrix(i,j)=1;
        end
    end
end
adMatrix=matrix;

end

