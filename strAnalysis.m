function [ T ] = strAnalysis( table )
fprintf('Decompose into Strong Connected Component!\n');
global iniTable;
%STRANALYSIS this function is to analysis the structure of the input table
%by using D-M decomposition
% It is the coarse decomposition , where structural
% Over-/Under-/Well-constrained subsystems are informed to the users;if the
% overconstraints are detected, we deleted them as default-line 54 & 55
%% find the unconnected subgraphs
H = height(table);
incmatrix = table2array(table);%extract table to matrix
adjmatrix = inc2adjm(incmatrix);%convert incidence matrix to adjacent matrix
% G=graph(adjmatrix);%not necessary
[n,s,m]=networkComponents(adjmatrix);%find the weakly connected components
fprintf('There are %d components in the system:\n',n);
%%
for i=1:n %extract each component from original table and generate the incidence matrix
    rows=[];
    fprintf('\n for component(%d)\n',i);
    for k=1:H
        if any(iniTable{k,m{i}}~= 0) 
            rows=[rows,k];
        end
    end
    t=iniTable(rows,m{i});
%     for ii=1:length(rows)
%         for jj=1:length(m{i})
%             if t{ii,jj}~= 0
%                 t{ii,jj}=1;
%             end
%         end
%     end
k=1;
%% for future usage of detecting redundant/conflicting group, here create a directed graph between varaibles and equations: variable--->equation
incmatrix = table2array(t);%extract table to matrix
[mm,nn]=size(incmatrix);
for i=1:mm
    for j=1:nn
        if incmatrix(i,j)~= 0
            incmatrix(i,j)=1;
        end
    end
end
V_V=zeros(nn,nn);V_E=zeros(nn,mm);
E_V=incmatrix;E_E=zeros(mm,mm);
Temp1=horzcat(V_V,V_E);
Temp2=horzcat(E_V,E_E);
Adjmatrix=vertcat(Temp1,Temp2);       
Adjmatrix=sparse(Adjmatrix);
names=[t.Properties.VariableNames,t.Properties.RowNames'];
compDAG=digraph(Adjmatrix,names);
tic
while ~isempty(t)
    fprintf('Structure Analysis Starts %dth time:\n',k);
%     T=timeit(dmAnalysis(t));
    [t,compDAG]=dmAnalysis(t,compDAG);
%     time=time+T;
%     fprintf('The running time is %d seconds\n',T);
%     k=k+1;
end
toc
% fprintf('Structure Analysis Starts %dth time:\n',k);
end
end
  




