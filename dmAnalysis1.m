function [finTable] = dmAnalysis1( T )
%DMANALYSIS structural decompose the input incidence matrix.
%   this function is based on dmperm function. Detail explanations are added.
% global iniTable;
% finTable=iniTable;
e=1e-11;
global iniTable;
finTable=iniTable(T.Properties.RowNames,T.Properties.VariableNames);
Over = table; A11=table;A12=table;
Under = table; A23=table;A34=table;
Well = table; A44=table;Block=table;
[m,n]=size(T);
for i=1:m
    for j=1:n
        if T{i,j}~= 0
            T{i,j}=1;
        end
    end
end
incmatrix=table2array(T);
[p,q,r,s,cc,rr] = dmperm(incmatrix);%decompose the adjacent matrix using DM algorithm
T=T(p,q);%Reorder the adjacent table 
%% structural analysis of the adjacent table after D-M decomposition
% coarse decomposition for over-/under-/well-constrained subparts-for
% details , refer to dmperm function of matlab
if rr(2)-1>= rr(1) && cc(2)-1>=cc(1)
    A11=T(rr(1):rr(2)-1,cc(1):cc(2)-1);
else
    fprintf('There is no unmatched variables!\n')
end
if   rr(1)<=rr(2)-1 && cc(2)<=cc(3)-1
    A12=T(rr(1):rr(2)-1,cc(2):cc(3)-1);
end
Under=[A11,A12];
if rr(3)<=rr(4)-1 && cc(4)<=cc(5)-1
    A34=T(rr(3):rr(4)-1,cc(4):cc(5)-1);
end
if rr(4)<=rr(5)-1 && cc(4)<=cc(5)-1
    A44=T(rr(4):rr(5)-1,cc(4):cc(5)-1);
end
Over=[A34;A44];
if rr(2)<=rr(3)-1 && cc(3)<=cc(4)-1
    A23=T(rr(2):rr(3)-1,cc(3):cc(4)-1);
end
Well=A23;
if isempty(Over)
    fprintf('There is no structural Over-constrained subsystem!\n');
else %if there exist overconstrained parts,analysis the structure by finding strongly connected components 
    fprintf('There exists structural Over-constrained subsystem!\n');
    [k,l]=size(A44);
    fprintf('The number is %d,they are: \n',k);
    for i=1:k
        fprintf(' %s', A44.Properties.RowNames{i})
    end
    %% Delete the structural overconstraints
    finTable(A44.Properties.RowNames,:)=[];
    fprintf('\n These overconstraints have been removed automaticly\n');
    
    %% Handle the strong connected components in the remaining part A34
    incmatrix=table2array(A34);
    [val,mi,mj]=bipartite_matching(incmatrix);% mi-->rows;mj-->columns
    [S,digraph,finTable]=vecomp(A34,mi,mj,finTable,'true');
    return
end
if isempty(Well)
    fprintf('There is no structural Well-constrained subsystem!\n');
else
    fprintf('There exists structural Well-constrained subsystem!\n');
    %%find the strong connected components along the hierarchical solving
    %%order
    incmatrix=table2array(A23);
    [val,mi,mj]=bipartite_matching(incmatrix);% mi-->rows;mj-->columns
    
    [S,digraph,finTable]=vecomp(A23,mi,mj,finTable,'true');
    return
    %% Or for extract the strong connected component in Well-constrained part, the following can be used(details infer to dmperm function in MATLAB)
%     fprintf('There are %d strong connected components(blocks) in this subsystem\n',S);%fine decomposition
%     if isempty(Under)
%         k=0;
%     else
%         k=1;
%     end
%      if isempty(Over)
%         kk=0;
%     else
%         kk=1;
%      end
%      j=1;
%     for i=1+k:length(r)-kk-1
%         fprintf('For block %d: \n',j);
%         Block = T(r(i):r(i+1)-1,s(i):s(i+1)-1);
%         Block = finTable(Block.Properties.RowNames,Block.Properties.VariableNames);
% %         Block=[Block,iniTable(Block.Properties.RowNames,end)];
%         [Block_af,overconstraints]=tablerref(Block,e);
%         if isempty(overconstraints)
%             fprintf('There is no numerical overconstraints in this Block!\n');%求解方程并且输入到以下Block

%         else
%             l=length(overconstraints);
%             fprintf('There are %d numerical overconstraints: \n',l);
%             for i=1:l
%                 fprintf(' %s', overconstraints{i});
%             end
%             finTable(overconstraints,:)=[];
%             fprintf('\nThese overconstraints have been removed automaticly\n');
%             return
%             %         finTable=Under(:,1:end-1);
%         end
%         j=j+1;
%     end
end
if isempty(Under) 
    fprintf('There is no structural Under-constrained subsystem!\n');
else 
    fprintf('There exists structural Under-constrained subsystem!\n');
%     Under=finTable(Under.Properties.RowNames,Under.Properties.VariableNames);%concatenate the last columns to the table
%     Under=[Under,iniTable(Under.Properties.RowNames,end)];
    incmatrix=table2array(Under);
    [val,mi,mj]=bipartite_matching(incmatrix);% mi-->rows;mj-->columns
    [S,C,finTable]=vecomp(Under,mi,mj,finTable,'true');
%     fprintf('There are %d strong connected components(blocks) in this subsystem\n',S);%fine decomposition
    Under=iniTable(Under.Properties.RowNames,Under.Properties.VariableNames);
    B=iniTable(Under.Properties.RowNames,end);
    [conflicting,redundant]=tablerref([Under,B]);
    if isempty(conflicting) && isempty(redundant)
        fprintf('There is no numerical overconstraints in this under-constrained subsystem!\n');
         finTable=[];%terminate the loop
    else
        if  ~isempty(conflicting)
            fprintf('\nThere exist conflicting constraints in Under-constrained system\n');
            l=length(conflicting);
            fprintf('The number is %d,they are: \n',l);
            for i=1:l
                fprintf(' %s ', conflicting{i});
            end
            finTable(conflicting,:)=[];
            fprintf('\nThese conflicting constraints have been removed automaticly\n');
        end
        if  ~isempty(redundant)
            fprintf('\nThere exist redundant constraints in Under-constrained system\n');
            l=length(redundant);
            fprintf('The number is %d,they are: \n',l);
            for i=1:l
                fprintf(' %s ', redundant{i});
            end
            finTable(redundant,:)=[];
            fprintf('\nThese redundant constraints have been removed automaticly\n');
        end
                
%         l=length(overconstraints);
%         fprintf('There are %d numerical overconstraints: \n',l);
%         for i=1:l
%             fprintf(' %s', overconstraints{i});
%         end
%         finTable(overconstraints,:)=[];
%         fprintf('\nThese overconstraints have been removed automaticly\n');
    end
   
end



end

