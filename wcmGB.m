function [ conflicting,redundant,xsolution] = wcmGB( scctPara,scctVar,compDAG,tol)
%wcmGB a combination of witness configuration method and grobner base to
%detect the conflicting and redundant constraints for non-linear system ; 
% witness configuration method are used to identify the overconstraints
% while grobner base are used further to distiguish between redundant and
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
%% apply WCM to detect the overconstraints
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
       overGroup=scctPara(r(i+1)+1:n,:);% over-constrained group of equations:table form:these constraints are redundant constraints
       over=vpa(table2array(overGroup));% over-constrained group of equations:matrix form:these constraints are redundant constriants
        break;
    end
end
[L_R,W_R]=size(over);% record the length of redundant constraints
%% check if there exist conflicting constriants in the basis group
chx={};
for j=1:length(varPara)
    chx=[chx,char(varPara(j))];
end
[mm,nn]=size(basis);
chbasis={};
for j=1:mm
    chbasis=[chbasis,char(basis(j))];
end
Rgd=groebner(chbasis,'lex',chx,tol);
while strcmp(Rgd,'1')
    temp=basis(end);
    tempT=basisGroup(end,:);
    basis(end)=[];% reset basis 
    basisGroup(end,:)=[];% reset basisgroup
    varPara=symvar(basis);
    chx={};
    for j=1:length(varPara)
        chx=[chx,char(varPara(j))];
    end
    [mm,nn]=size(basis);
    chbasis={};
    for j=1:mm
        chbasis=[chbasis,char(basis(j))];
    end
    Rgd=groebner(chbasis,'lex',chx,tol);
    over=[over;temp];% add this equation to over:this is conflicting equations
    overGroup=[overGroup;tempT];% add this equation to overgroup:this equation is conflicting
    if  ~strcmp(Rgd,'1')
        break;
    end
end
%% testing conflicting and redundant constraints
xsolution=x(locb);
if  isempty(over)            
    return;
else
    [MM,nn]=size(over);
    for i=1:MM
%         varOver=symvar(over(i));% find the varaibles of overconstraints
%         %     varOver=symvar(over);
%         eNames=[];
%         for j=1:length(varOver)
%             eNames=[eNames,predecessors(compDAG,char(varOver(j)))'];
%         end
%         eNames=unique(eNames);
%         eNames=intersect(eNames,scctPara.Properties.RowNames);
%         eNames= setdiff(eNames,overGroup(i,:).Properties.RowNames);
%         basisGroup=scctPara(eNames,:);% redefine the basis group that are conflicting/redundant with over(i);
%         basis=vpa(table2array(basisGroup));% redefine the basis group that are conflicting/redundant with over(i);
        group=basisGroup.Properties.RowNames;% rownames of basis group
        source_group=findSource(compDAG,group);
        %% apply GB to distinguish redundant and conflicting constraints
%         [mm,nn]=size(basis);
%         chbasis={};
%         for j=1:mm
%             chbasis=[chbasis,char(basis(j))];
%         end
%         chx={};
%         for j=1:length(varPara)
%             chx=[chx,char(varPara(j))];
%         end
%         Rgd=groebner(chbasis,'lex',chx,tol);
%         % formatSpec = 'groebner::gbasis([%s], LexOrder)';
%         % str = sprintf('%s,',chbasis{:});
%         % str(end)=[];
%         % str = sprintf(formatSpec,str);
%         % Rgd=evalin(symengine,str);
%         %     [MM,nn]=size(over);
%         %     for i=1:MM
%         testOver=char(over(i)); % Drop the idea of using Grobner bases
%         chbasis=[chbasis,testOver];% to test conflicting/redundant
%         rgd=groebner(chbasis,'lex',chx,tol);% apply grobner basis to find reduced GB
        %             formatSpec = 'groebner::gbasis([%s], LexOrder)';
        %             str = sprintf('%s,',chbasis{:});
        %             str(end)=[];
        %             str = sprintf(formatSpec,str);
        %             rgd=evalin(symengine,str);
        if i > L_R % replace strcmp(rgd,'1')
            cName=overGroup(i,:).Properties.RowNames;
            source_C = findSource(compDAG,cName);
            Group=unique([group',source_group,source_C]);
            if ~isempty(Group)
                disp([cName,'is conflicting with group:',Group]);
            end
            conflicting=[conflicting,cName];
        else                                  %if  isempty(setdiff(rgd,Rgd))
            rName=overGroup(i,:).Properties.RowNames;
            source_C = findSource(compDAG,rName);
            Group=unique([group',source_group,source_C]);
            if ~isempty(Group)
                disp([rName,'is redundant with group:',Group]);
            end
            redundant=[redundant,rName];
%         else
%             error('Error. \n This constraint is neither conflicting nor redundant but it is detected by WCM as a redundant/conflicting constraint!')
        end
    end
end
end
%%
    




