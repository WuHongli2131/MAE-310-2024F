function [str,stxi,tor,xita]=stresspoly(T,R,x,y)%把
xita=atan2(y,x);
r=sqrt(x^2+y^2);
str=T/2*(1-R^2/r^2)+T/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*xita);
stxi=T/2*(1+R^2/r^2)-T/2*(1+3*R^4/r^4)*cos(2*xita);
tor=-T/2*(1+2*R^2/r^2-3*R^4/r^4)*sin(2*xita);
end
%验证完毕，没有问题