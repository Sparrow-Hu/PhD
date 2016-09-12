function plotParapts( nurbsurface,parametricpoints )
%PLOTPARAPTS plot the parametric points on a nurbsurface
%   nurbsurface is a nurbs surface; parametric points are U,V parameters in
%   parametric space
srf=nurbsurface;
[m,n]=size(parametricpoints);
for i=1:m
    temp=parametricpoints(i,:);
    pp=nrbeval(srf,{temp(1),temp(2)});
    X(i)=pp(1);
    Y(i)=pp(2);
    Z(i)=pp(3);
end
plot3(X,Y,Z,'b+','MarkerSize',20);

end

