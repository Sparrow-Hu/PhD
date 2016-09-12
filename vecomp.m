function [conflicting,redundant,finTable,compDAG] = vecomp( table,finTable,compDAG )
%VECOMP Handle the strong connected components in table. If it is the full
%rank, then it is solved. If there exist overconstraints, then the
%overconstraints are deleted! Update the corresponding finTable. finTable
%is adjacent table;table is the subpart;fintable is the coefficient
%version of finTable
global iniTable;
tol=1e-4;conflicting=[];redundant=[];
fintable=iniTable(finTable.Properties.RowNames,finTable.Properties.VariableNames);%transform incident table to 
iniTable{fintable.Properties.RowNames,end}=-iniTable{fintable.Properties.RowNames,end};
B=iniTable(fintable.Properties.RowNames,end);
%% Apply direct numerical method on this subpart
% tic
% fprintf('Direct Detection using QR factorization: \n');
% [Conflicting,Redundant]=qrRC([fintable,B],compDAG,tol)
% toc
% fprintf('Direct Detection using LU factorization: \n');
% [Conflicting,Redundant]=luRC([fintable,B],1e-4)
% tic
% fprintf('Direct Detection using Gauss-Jordan factorization: \n');
% [Conflicting,Redundant]=tablerref([fintable,B],1e-4)
% toc
while ~isempty(table)
    [dgSC,s]=dagScomp(table);% dag(dgSC) of strong connected components(s)
    %% Identify the overconstraints in SCC which does not have predecessors
%     plot(dgSC);
    l=length(s);k=0;
    for i=1:l
        preIDs=predecessors(dgSC,i)';
        if isempty(predecessors(dgSC,i)) && length(s{i})>1
            temp=s{i};
            k=k+1;% a flag to terminate the loop
            L=length(temp);
            sccT=iniTable(temp(L/2+1:L),temp(1:L/2));
            B=iniTable(sccT.Properties.RowNames,end);
            %         [conflicting,redundant] = tablerref([sccT,B],tol);%First trying:Gauss-Jordan Elimination with rref:However it is not stable
            %         num_os=tableRC([sccT,B]);%Second trying: Comparing SVD and LU; LU is better but still bad
            [conflicting,redundant]=qrRC([sccT,B],compDAG,tol);%Thrid trying:SVD&QR&LU, we do not specify the tolerence but classify distance
            if isempty(conflicting) && isempty(redundant)
                fprintf('There is no numerical overconstraints in Block %d\n',i);
                %% based on compDAG, we add some edges from equations to variables in this block: for each variable, we add edges from all the equations to this variable
                [mm,nn]=size(sccT);
                for i=1:nn
                    for j=1:mm
                        compDAG=addedge(compDAG,sccT.Properties.VariableNames{i},sccT.Properties.RowNames{j},1);% so that we can know, for a solution X, it is solved from which equations
                    end
                end
                coeMatrix=table2array(sccT);
                b=table2array(B);
                X=linsolve(coeMatrix,b);%get the solution
                %% input the solution and update the finTable
                [m,n]=size(finTable);
                B=iniTable{finTable.Properties.RowNames,end};%reset the B for new use
                [lia,locb]=ismember(finTable.Properties.VariableNames,sccT.Properties.VariableNames);
                for i=1:n%columns
                    if lia(i)
                        for j=1:m %update the b value of this subpart(without structural overconstraints)
                            B(j)=B(j)-fintable{j,i}*X(locb(i));
                        end
                    end
                end
                iniTable{finTable.Properties.RowNames,end}=B;
                finTable(temp(L/2+1:L),:)=[];
                finTable(:,temp(1:L/2))=[];
                table(temp(L/2+1:L),:)=[];
                table(:,temp(1:L/2))=[];
                fintable(temp(L/2+1:L),:)=[];
                fintable(:,temp(1:L/2))=[];
                fprintf('\nThis Block has been released\n');
            end
            if  ~isempty(conflicting)
                fprintf('There exist conflicting constraints in Block %d\n',i);
                l=length(conflicting);
                fprintf('The number is %d,they are: \n',l);
                for i=1:l
                    fprintf(' %s ', conflicting{i});
                end
                finTable(conflicting,:)=[];
                table(conflicting,:)=[];
                fintable(conflicting,:)=[];
                fprintf('These conflicting constraints are removed\n');
            end
            if  ~isempty(redundant)
                fprintf('There exist redundant constraints in Block %d\n',i);
                l=length(redundant);
                fprintf('The number is %d,they are: \n',l);
                for i=1:l
                    fprintf(' %s ', redundant{i});
                end
                finTable(redundant,:)=[];
                table(redundant,:)=[];
                fintable(redundant,:)=[];
                fprintf('These redundant constraints are removed\n');
            end
        end
    end
    [mm,nn]=size(finTable);
    r=[];c=[];
    for i=1:mm
        if all(finTable{i,:}==0)
            temp=finTable(i,:);
            d=iniTable{temp.Properties.RowNames,end};
            name=char(temp.Properties.RowNames);
            fprintf('After solving,the difference for %s is %d\n',name,d);
            group = findSource(compDAG,{name});% filter the variables and keep the equations in V
            if abs(d)<= tol
                fprintf('Thus %s is redundant\n',name); 
                r=[r,temp.Properties.RowNames];
                disp([name,'is redundant with group:',group]);
            else
                fprintf('Thus %s is conflicting\n',name); 
                c=[c,temp.Properties.RowNames];
                disp([name,'is conflicting with group:',group]);
            end
        end
    end
    redundant=[redundant,r];
    if ~isempty(r)
        finTable(r,:)=[];
        table(r,:)=[];
        fintable(r,:)=[];
        fprintf('The redundant constraints obtained by feeding process are removed\n');
    end
    conflicting=[conflicting,c];
    if ~isempty(c)
        finTable(c,:)=[];
        table(c,:)=[];
        fintable(c,:)=[];
        fprintf('The conflicting constraints obtained by feeding process are removed\n');
    end
    if k==0%if there is no frontier strong connected components(components that have no predecessors),exit the loop
        break;
    end
end

end

    

               
               
        




