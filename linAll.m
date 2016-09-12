function [ iniTable ] = linAll( ceqPara,ceqVars )
%linearize the nonlinear equations to first order
%   output Ax=b in table form
[m,~]=size(ceqPara);
if m ~= 1
    ceqPara=ceqPara.';% rows to columns
    ceqVars=ceqVars.';% rows to columns
end
varPara=symvar(ceqPara);% find the variables of ceqPara
varVars=symvar(ceqVars);% find the variables of ceqVars
f=symfun(ceqVars,varVars);
v=@(X)double(subs(f,varVars,X));
%% generate witness configurations
r(1)=1;
for i=1:10% find the number of basis constraints
    options=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',1000,'Display','iter','PrecondBandWidth',1);
    [x,fval,exitflag,output] = fsolve(v,yRandom(length(varVars),10,1)',options);% start with radomly permutated numbers with 1 as mean and stand deviation
    [Lia,locb]=ismember(varPara,varVars);
    Jamatrix=jacobian(ceqPara,varPara);
    Jamatrix=subs(Jamatrix,varPara,x(locb));
    r(i+1)=rank(Jamatrix);
    if r(i+1)==r(i)
%         Jamatrix=Jamatrix';
%         [~,n]=size(Jamatrix);
%         [~,~,E] = qr(Jamatrix,0);
%         ceqPara=ceqPara(E,:);   %adjust the order of constraints equations to pump  the overconstraints to the bottom
%         basisGroup=ceqPara(1:r(i+1),:);% basis group of equations:table form
%         basis=vpa(table2array(basisGroup));% basis group of equations:matrix form
%         overGroup=ceqPara(r(i+1)+1:n,:);% over-constrained group of equations:table form
%         over=vpa(table2array(overGroup));% over-constrained group of equations:matrix form
        break;
    end
end
%% linearize to first order using taylor series
xsolution=x(locb);
linPara=taylor(ceqPara, varPara, xsolution, 'Order', 2);
b= -double(subs(linPara,varPara, zeros(1,length(varPara))));
Jamatrix=double(jacobian(linPara,varPara));
varCell={};
for i=1:length(varPara)
    varCell=[varCell,char(varPara(i))];
end
varCell=[varCell,'b'];
iniTable =array2table([Jamatrix,b'],'VariableNames',varCell); 

end

