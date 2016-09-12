function [T,compDAG] = dmAnalysis( T,compDAG)
%DMANALYSIS structural decompose the input incidence matrix.
%   this function is based on dmperm function. Detail explanations are added.
% global iniTable;
% finTable=iniTable;
tol=1e-4;
global iniTable ceqPara ceqVars;
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
    fprintf('\n Handling the over-constrained part!\n');
    if isempty(ceqPara)
        [~,~,T,compDAG]=vecomp(Over,T,compDAG);% linear case
    else
        [~,~,T,compDAG]=vecompn(Over,T,compDAG);% nonlinear case
    end
    return
end
if isempty(Well)
    fprintf('There is no structural Well-constrained subsystem!\n');
else
    fprintf('There exists structural Well-constrained subsystem!\n');
    fprintf('\n Handling the Well-constrained part!\n');
    if isempty(ceqPara)
        [conflicting,redundant,T,compDAG]=vecomp(Well,T,compDAG);% linear case
    else
        [conflicting,redundant,T,compDAG]=vecompn(Well,T,compDAG);% nonlinear case
    end
    if isempty(Under) && isempty(conflicting) && isempty(redundant)
        T=[];  
    end
    return
end
if isempty(Under)
    fprintf('There is no structural Under-constrained subsystem!\n');
else
    fprintf('There exists structural Under-constrained subsystem!\n');
    fprintf('\n Handling the Under-constrained part!\n');
    if isempty(ceqPara)
        [Conflicting,Redundant,T,compDAG]=vecomp(Under,T,compDAG);% linear case
        fprintf('G-E on the remaining part!\n');
        finTable=iniTable(T.Properties.RowNames,T.Properties.VariableNames);
        B=iniTable(T.Properties.RowNames,end);
        [conflicting,redundant]=qrRC([finTable,B],compDAG,tol);%Second trying:QR factorization
    else
        [Conflicting,Redundant,T,compDAG]=vecompn(Under,T,compDAG);% non-linear case
        fprintf('WCM on the remaining part!\n');
        comp_ceqPara=ceqPara(T.Properties.RowNames,:);
        comp_ceqVars=ceqVars(T.Properties.RowNames,:);
        [m,n]=size(comp_ceqPara);
        if m<10% if there are less than 10 equations,Use GB
             fprintf('Detection using Grobner Basis!');
            [conflicting,redundant]=wcmGB(comp_ceqPara,comp_ceqVars,compDAG,tol);%apply wcm and GB to find redundant and conflicting constraints
        else% else use OPtimization
            fprintf('Detection using Optimization!');
            [conflicting,redundant]=wcmOP(comp_ceqPara,comp_ceqVars,compDAG,tol);%apply wcm and OP to find redundant and conflicting constraints
        end
    end
    %% output results
    if  ~isempty(conflicting)
        fprintf('There exist conflicting constraints in Block %d\n',i);
        l=length(conflicting);
        fprintf('The number is %d,they are: \n',l);
        for i=1:l
            fprintf(' %s ', conflicting{i});
        end
        T(conflicting,:)=[];
        fprintf('These conflicting constraints are removed\n');
    end
    if  ~isempty(redundant)
        fprintf('There exist redundant constraints in Block %d\n',i);
        l=length(redundant);
        fprintf('The number is %d,they are: \n',l);
        for i=1:l
            fprintf(' %s ', redundant{i});
        end
        T(redundant,:)=[];
        fprintf('These redundant constraints are removed\n');
    end
    if isempty(conflicting) && isempty(redundant) && isempty(Conflicting) && isempty(Redundant)
        fprintf('There is no numerical overconstraints in this under-constrained subsystem!\n');
        T=[];%terminate the loop
    end   
end
end

