function [ xx,fval,varscstr] = OP( basis,over,x,vars )
%OP Summary of this function goes here
%   Detailed explanation goes here
%     function [c,ceq]=nonlinear(X)
%         c=[];
%         ff=symfun(basis,varscstr);% it is a temporary function for the usage of fsolve and it is necessary and important
%         ceq=double(subs(ff,varscstr,X));
%     end
% define objective function
% obj=sumsqr([basis,over]);
obj=0;
for i=1:length(basis)
    obj=obj+basis(i)^2;
end
obj=obj+over^2;
varscstr = symvar(obj);% find the symbolic variables in these symbolic expressions
F=symfun(obj,varscstr);
[Lia,locb]=ismember(varscstr,vars);
objective=@(X)double(subs(F,varscstr,X));
% initialize parameters for optimization
% A=[];b=[];Aeq=[];beq=[];lb=[];ub=[];
% oldoptions=optimoptions(@fmincon,'MaxIter',1000,'MaxFunEvals',1000,'Display','iter');
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
for i=1:length(Lia)
    if Lia(i)==0
        x0(i)=0;
    else
        x0(i)=x(locb(i));
    end
end
problem = createOptimProblem('fmincon','x0',x0,'objective',objective,'options',opts);
gs = GlobalSearch;
[xx,fval] = run(gs,problem);
% [xx,fval,exitflag]=fmincon(objective,x0,A,b,Aeq,beq,lb,ub,@nonlinear,oldoptions);
% k=0;
% while exitflag==0
%     oldoptions.MaxIter=4000;
%     oldoptions.MaxFunEvals=4000;
%     [xx,fval,exitflag]=fmincon(objective,yRandom(length(x0),10,1)',A,b,Aeq,beq,lb,ub,@nonlinear,oldoptions);
%     k=k+1;
%     if k>5
%         break;
%     end
% end
% if exitflag==0
%     oldoptions.MaxIter=20000;
%     oldoptions.MaxFunEvals=40000;
%     [xx,fval,exitflag]=fmincon(objective,x(locb),A,b,Aeq,beq,lb,ub,@nonlinear,oldoptions);
% end
end

