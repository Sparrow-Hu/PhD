function [ conflicting,redundant,xsolution] = wcmOP( scctPara,scctVar,compDAG,tol)
%wcmGB a combination of witness configuration method and optimization method to
%detect the conflicting and redundant constraints for non-linear system ;
% witness configuration method are used to identify the overconstraints
% while optimization method are used further to distiguish between redundant and
% conflicting constraints.
% scctPara is the original constraint system where parameters are kept;
% scctVars is the modified version system where parameters are adapted to
% variables so that WCM method can be used
conflicting={};
redundant={};
paraArray=table2array(scctPara);
varPara=symvar(paraArray);% find the variables of paraArray
varArray=table2array(scctVar);
varVars=symvar(varArray);% find the variables of varArray
%% apply WCM to detect the overconstraints--that is , the redundant constraints
%     function v=V(X)% this function is created to adopt X as input since fsolve only accept one input;have to assemble the variables into this vector
%         f=symfun(varArray,varVars);% it is a temporary function for the usage of fsolve and it is necessary and important
%         v=double(subs(f,varVars,X));
%     end
f=symfun(varArray,varVars);
v=@(X)double(subs(f,varVars,X));% this is a kind of hack that can be used to replace matlabFunction for seperate outputs in case of input with one vector containing several variables
r(1)=1;
for i=1:10% find the number of basis constraints
    options=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',1000,'Display','iter','PrecondBandWidth',1);
    [x,fval,exitflag,output] = fsolve(v,yRandom(length(varVars),10,1)',options);% start with radomly permutated numbers with 1 as mean and stand deviation
    [Lia,locb]=ismember(varPara,varVars);
    Jamatrix=jacobian(paraArray,varPara);
    Jamatrix=subs(Jamatrix,varPara,x(locb));
    r(i+1)=rank(Jamatrix);
    if r(i+1)==r(i)
        Jamatrix=Jamatrix';
        [~,n]=size(Jamatrix);
        [~,~,E] = qr(Jamatrix,0);
        scctPara=scctPara(E,:);   %adjust the order of constraints equations to pump  the overconstraints to the bottom
        basisGroup=scctPara(1:r(i+1),:);% basis group of equations:table form
        basis=vpa(table2array(basisGroup));% basis group of equations:matrix form
        overGroup=scctPara(r(i+1)+1:n,:);% over-constrained group of equations:table form
        over=vpa(table2array(overGroup));% over-constrained group of equations:matrix form
        break;
    end
end
%% check if there exist conflicting constriants in the basis group// drop  the idea, since it is very time-consuming when the number of equations is bigger than 10
% chx={};
% for j=1:length(varPara)
%     chx=[chx,char(varPara(j))];
% end
% [mm,nn]=size(basis);
% chbasis={};
% for j=1:mm
%     chbasis=[chbasis,char(basis(j))];
% end
% Rgd=groebner(chbasis,'lex',chx,tol);
% while strcmp(Rgd,'1')
%     temp=basis(end);
%     tempT=basisGroup(end,:);
%     basis(end)=[];% reset basis and basis group
%     basisGroup(end,:)=[];
%     varPara=symvar(basis);
%     chx={};
%     for j=1:length(varPara)
%         chx=[chx,char(varPara(j))];
%     end
%     [mm,nn]=size(basis);
%     chbasis={};
%     for j=1:mm
%         chbasis=[chbasis,char(basis(j))];
%     end
%     Rgd=groebner(chbasis,'lex',chx,tol);
%     if  ~strcmp(Rgd,'1')
%         over=[over;temp];% add this equation to basis and basisgroup
%         overGroup=[overGroup;tempT];
%         break;
%     end
% end
%% testing conflcting constraints using optimization method: check if some of the redundant constraints are redudant, if yes, this means that there is no conflicting constraints else there exists conflicting constraints in the basis
xsolution=x(locb);
if  ~isempty(over)
    [MM,nn]=size(over);
    xx=x(locb);
    var=varPara;
    for i=1:MM
        %     varOver=symvar(over(i));% find the varaibles of overconstraints
        %     %     varOver=symvar(over);
        %     eNames=[];
        %     for j=1:length(varOver)
        %         eNames=[eNames,predecessors(compDAG,char(varOver(j)))'];
        %     end
        %     eNames=unique(eNames);
        %     eNames=intersect(eNames,scctPara.Properties.RowNames);
        %     eNames= setdiff(eNames,overGroup(i,:).Properties.RowNames);
        %     basisGroup=scctPara(eNames,:);% redefine the basis group that are conflicting/redundant with over(i);
        %     basis=vpa(table2array(basisGroup));% redefine the basis group that are conflicting/redundant with over(i);
        group=basisGroup.Properties.RowNames;% rownames of basis group
        source_group=findSource(compDAG,group);
        %% apply Optimization to distinguish redundant and conflicting constraints
        % define constraints---formed by basis constraints ;Here we have to switch
        % nested function for the output of non-linear constraints
        % substitue the values to overGroup and distinguish the redundant and
        % conflicting constraints
        % varscstr = symvar(basis);% find the symbolic variables in these symbolic expressions
        %     function [c,ceq]=nonlinear(X)
        %         c=[];
        %         ff=symfun(basis,varscstr);% it is a temporary function for the usage of fsolve and it is necessary and important
        %         ceq=double(subs(ff,varscstr,X));
        %     end
        % % define objective function
        % [MM,nn]=size(over);
        % obj=0;
        % for i=1:MM
        %     obj=obj+abs(over(i))^2;
        % end
        % F=symfun(obj,varscstr);
        % % objective=@(X)double(subs(F,varscstr,X));
        % objective=@(x)0;
        % % initialize parameters for optimization
        % A=[];b=[];Aeq=[];beq=[];lb=[];ub=[];
        % [Lia,locb]=ismember(varscstr,varVars);
        % oldoptions=optimoptions(@fmincon,'MaxIter',100000,'MaxFunEvals',100000,'Display','iter');
        % [xx,fval,exitflag]=fmincon(objective,x(locb),A,b,Aeq,beq,lb,ub,@nonlinear,oldoptions);
        % xsolution=xx;
        % if exitflag==0
        %     oldoptions.MaxIter=20000;
        %     oldoptions.MaxFunEvals=40000;
        %     [xx,fval,exitflag]=fmincon(objective,x(locb),A,b,Aeq,beq,lb,ub,@nonlinear,oldoptions);
        % end
        [ xx,fval,var] = OP1( basis,over(i),xx,var);
        %         if exitflag==0% optimization failed,possibly due to the existance of conflicting constraints in the basis constraints. Check again by other redundant constraints
        %             continue
        %         end
        if fval> tol
            %             error('Error. WCM detect conflicting constraints!');
            cName=overGroup(i,:).Properties.RowNames;
            source_C = findSource(compDAG,cName);
            Group=unique([group',source_group,source_C]);
            if ~isempty(Group)
                disp([cName,'is conflicting with group:',Group]);
            end
            conflicting=[conflicting,cName];
        else
            rName=overGroup(i,:).Properties.RowNames;
            source_C = findSource(compDAG,rName);
            Group=unique([group',source_group,source_C]);
            if ~isempty(Group)
                disp([rName,'is redundant with group:',Group]);
            end
            redundant=[redundant,rName];
        end
    end
    %     if ~isempty(conflicting)
    %         error('Error. WCM detect conflicting constraints!')
    %     end
%     if isempty(redundant)% if WCM detect overconstraints(redundant constraints) and we recheck them without finding any redundant constraints, meaning that there are conficting constraints in the basis ,resulting basis constraints unsolvable
%         fprintf('There exists conflicting constraints in the basis constraints!');
%         % check the constraint from the bottom one by one
%         [L,W]=size(basis);
%          j=1;
%         for i=1:L% this is the very safe way to check the conflicting constriants! if in the future , no more errors happen, using just for loop is enough!
%             [ xx,fval,var,exitflag] = OP( basis(1:L-i),basis(end-i+1),xx,var);
%             if exitflag==0 % optimnization failed ,there are more conflicting constraints!!!
%                 Cflicts2check(j)=basis(end-i+1);
%                 j=j+1;
%                 continue % there are more conflicting constraints in the basis(1:L-i) constraints! Check in the next loop!
%             end
%             Basis=basis(1:L-i);% find Basis constraints that independent
%             Cflicts2check(j)=basis(end-i+1);% find potential conflicting constraints to be checked
%             break
%         end
%         if ~isempty(Cflicts2check)
%             for i=1:j
%                 [ xx,fval,var,exitflag] = OP( Basis,Cflicts2check(i),xx,var);
%                 if abs(fval)> tol
%                     cName=overGroup(i,:).Properties.RowNames;
%                     source_C = findSource(compDAG,cName);
%                     Group=unique([group',source_group,source_C]);
%                     if ~isempty(Group)
%                         disp([cName,'is conflicting with group:',Group]);
%                     end
%                     conflicting=[conflicting,cName];
%                 else
%                     error('Error. Optimization method find redundant constraints in the basis!');
%                     %                     rName=overGroup(i,:).Properties.RowNames;
%                     %                     source_C = findSource(compDAG,rName);
%                     %                     Group=unique([group',source_group,source_C]);
%                     %                     if ~isempty(Group)
%                     %                         disp([rName,'is redundant with group:',Group]);
%                     %                     end
%                     %                     redundant=[redundant,rName];
%                 end
%                 %                 if ~isempty(redundant)
%                 %                 end
%             end
%             
%         end
%     end
end


