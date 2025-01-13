function [x,y]=normvector(u,v)
x=-v/sqrt(u^2+v^2);
y=u/sqrt(u^2+v^2);
end
%判断法向方向需要在边界上处理