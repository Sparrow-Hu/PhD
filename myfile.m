function out1 = myfile(x,y,z)
%MYFILE
%    OUT1 = MYFILE(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    13-Jun-2016 23:26:52

t2 = x.^2;
t3 = y.^2;
t4 = z.^2;
t5 = t2+t3+t4;
out1 = log(t5)+1.0./sqrt(t5);