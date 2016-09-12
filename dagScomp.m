function [ dgSC,s ] = dagScomp( table )
%DAGSCOMP is to draw DAG between strong connnected components in table
%  table is a adjacent table;dgSC is the DAG of s while s is a cell
%  containing strong connected components
incmatrix=table2array(table);
[val,mi,mj]=bipartite_matching(incmatrix);% mi-->rows;mj-->columns;maximum matching between incidence matrix
[m,n]=size(incmatrix);
V_E=incmatrix;E_V=zeros(n,m);
for i=1:length(mj)
   E_V(mj(i),mi(i))=1;
end
V_V=zeros(n,n);E_E=zeros(m,m);
Temp1=horzcat(V_V,E_V);
Temp2=horzcat(V_E,E_E);
Adjmatrix=vertcat(Temp1,Temp2);       
Adjmatrix=sparse(Adjmatrix);
names=[table.Properties.VariableNames,table.Properties.RowNames'];
dg=digraph(Adjmatrix,names);
s = conncomp(dg,'OutputForm','cell');
l=length(s);% construct dag of strong connected components
adj=zeros(l,l);
for i=1:l
    preIDs={};% the predecessors of this component s{i}
    ss=s{i};
    for k=1:length(ss)
        preIDs=[preIDs;predecessors(dg,ss(k))];
    end
    for j=1:l
        if i ~= j
            z=ismember(preIDs,s{j});
            if any(z)
                adj(i,j)=1;%in the future if need to know how many edges are involved , just summerzie the 1s in z; current, we do not consider num(edges) 
            end
        end
    end
end
dgSC=digraph(adj);
% plot(dgSC)

end

