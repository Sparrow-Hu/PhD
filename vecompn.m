function [conflicting,redundant,finTable,compDAG] = vecompn( table,finTable,compDAG )
%VECOMP Handle the strong connected components in table. If it is the full
%rank, then it is solved. If there exist overconstraints, then the
%overconstraints are deleted! Update the corresponding finTable. finTable
%is adjacent table;table is the subpart;fintable is the coefficient
%version of finTable
global iniTable ceqPara ceqVars ;
tol=1e-4;conflicting=[];redundant=[];
comp_ceqPara=ceqPara(finTable.Properties.RowNames,:);% equations with respect to equations in the component(parameters are fixed)
comp_ceqVars=ceqVars(finTable.Properties.RowNames,:);% equations with respect to equations in the component(parameters are free)
%% Apply direct numerical method on this subpart
% tic
% fprintf('Direct Detection using wcm and OP: \n');
% [Conflicting,Redundant]=wcmOP(comp_ceqPara,comp_ceqVars,compDAG,tol);
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
            scctVar=ceqVars(temp(L/2+1:L),:);
            scctPara=ceqPara(temp(L/2+1:L),:);
            [m,n]=size(scctPara);
            if m<10% if there are less than 10 equations,Use GB
                fprintf('Detection using Grobner Basis!');
                [conflicting,redundant,xsolution]=wcmGB(scctPara,scctVar,compDAG,tol);%apply wcm and GB to find redundant and conflicting constraints
            else% else use OPtimization
                fprintf('Detection using Optimization!');
                [conflicting,redundant,xsolution]=wcmOP(scctPara,scctVar,compDAG,tol);%apply wcm and OP to find redundant and conflicting constraints
            end
            if isempty(conflicting) && isempty(redundant)
                fprintf('There is no numerical overconstraints in Block %d\n',i);
                %% based on compDAG, we add some edges from equations to variables in this block: for each variable, we add edges from all the equations to this variable
                [mm,nn]=size(sccT);
                for i=1:nn
                    for j=1:mm
                        compDAG=addedge(compDAG,sccT.Properties.VariableNames{i},sccT.Properties.RowNames{j},1);% so that we can know, for a solution X, it is solved from which equations
                    end
                end
                paraArray = table2array(scctPara);
                vars = symvar(paraArray);% find the symbolic variables in these symbolic expressions
                f=symfun(paraArray,vars);% it is a temporary function for the usage of fsolve and it is necessary and important
                v=@(X)double(subs(f,vars,X));
                %                 for i=1:3
                options = optimoptions(@fsolve,'MaxIter',1000,'MaxFunEvals',1000,'Display','iter','PrecondBandWidth',1);
                [X,exitflag] =fsolve(v,xsolution,options);% get the solution
                %                     if exitflag==-1 | exitflag==-2 | exitflag==-3
                %                         options = optimoptions(@fsolve,'Algorithm','trust-region-reflective','MaxIter',1000,'MaxFunEvals',1000,'Display','iter','PrecondBandWidth',1);
                %                         [X,fval,exitflag,output] =fsolve(v,X,options);% get the solution
                %                     end
                %                     if exitflag==-1 | exitflag==-2 | exitflag==-3
                %                         options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',1000,'Display','iter','PrecondBandWidth',1);
                %                         [X,fval,exitflag,output] =fsolve(v,X,options);% get the solution
                %                     end
                while exitflag==0
                    options.MaxIter=10000;
                    options.MaxFunEvals=10000;
                    [X,exitflag]= fsolve(v,yRandom(length(vars),10,1)',options);
                end
                
                %% back substitution to the original equations and update adjacent table/ceqPara/ceqVars
                [m,n]=size(ceqPara);
                for i=1:m % substitute to the initial equations and update
                    tempPara=ceqPara{i,:};
                    tempVar=ceqVars{i,:};
                    tempPara=subs(tempPara,vars,X);
                    tempVar=subs(tempVar,vars,X);
                    ceqPara{i,:}=tempPara;
                    ceqVars{i,:}=tempVar;
                end
                finTable(temp(L/2+1:L),:)=[];
                finTable(:,temp(1:L/2))=[];
                table(temp(L/2+1:L),:)=[];
                table(:,temp(1:L/2))=[];
                ceqPara(temp(L/2+1:L),:)=[];
                ceqVars(temp(L/2+1:L),:)=[];
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
                ceqPara(conflicting,:)=[];
                ceqVars(conflicting,:)=[];
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
                ceqPara(redundant,:)=[];
                ceqVars(redundant,:)=[];
                fprintf('These redundant constraints are removed\n');
            end
        end
    end
    [mm,nn]=size(ceqPara);
    r=[];c=[];
    for i=1:mm
        temp=ceqPara(i,:);
        t=table2array(temp);
        if isempty(symvar(t))
            d=double(t);
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
        ceqPara(r,:)=[];
        ceqVars(r,:)=[];
        fprintf('The redundant constraints obtained by feeding process are removed\n');
    end
    conflicting=[conflicting,c];
    if ~isempty(c)
        finTable(c,:)=[];
        table(c,:)=[];
        ceqPara(c,:)=[];
        ceqVars(c,:)=[];
        fprintf('The conflicting constraints obtained by feeding process are removed\n');
    end
    if k==0%if there is no frontier strong connected components(components that have no predecessors),exit the loop
        break;
    end
end

end
