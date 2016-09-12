function [ group ] = findSource( digraph,enames)
%FINDSOURCE Given a set of equations, finding their source equations: if
%some of its variable has been solved, find the equations who give such
%solution
%   Digraph is a directed graph which initially all equations-->variables; if
%   some variables are solved by equations, then these variables-->equations
global iniTable
group={};
L = length(enames);
for i=1:L
    v=dfsearch(digraph,char(enames(i)));% trace back all the variables and equations using depth first search
    v(1)=[];
    temp = intersect(v,iniTable.Properties.RowNames);% filter the variables and keep the equations in V
    group=[group,temp'];
end


end

