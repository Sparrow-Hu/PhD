% syms x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36 x37 x38 x39 x40 x41


A=[];
b=[];
temp_Aeq=T(finTable.Properties.RowNames,finTable.Properties.VariableNames);
Aeq=table2array(temp_Aeq);
temp_beq=T(finTable.Properties.RowNames,end);
beq=table2array(temp_beq);
x0=linsolve(Aeq,-beq);
[x,fval] = fmincon(@Objective,x0,A,b,Aeq,-beq);