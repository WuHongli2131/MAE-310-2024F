%这一部分主要是为了极坐标转化成笛卡尔
function [stx sty tau]=polytocoor(str,stxi,tor,xita)
stx= (str+stxi)/2+(str-stxi)/2*cos(2*xita)+tor*sin(2*xita);
sty=(str+stxi)/2-(str-stxi)/2*cos(2*xita)-tor*sin(2*xita);
tau=(-str-stxi)/2*sin(2*xita)+tor*sin(2*xita);
end