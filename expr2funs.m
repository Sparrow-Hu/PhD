function [ fun ] = expr2funs( symArrays )
%EXPR2FUNS convert symbolic arrays to functions
%   by looping on matlab buit-in function matlabFunction
[m,n]=size(symArrays);
for i=1:m
    fun(i)=sym(matlabFunction(symArrays(i)));
end
end

