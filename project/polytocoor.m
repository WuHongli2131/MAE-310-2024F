%这一部分主[]要是为了极坐标转化成笛卡尔
function [stx, sty, tau]=polytocoor(str,stxi,tor,xita)
stx= str*cos(-xita)^2+stxi*sin(-xita)^2+2*tor*sin(-xita)*cos(-xita);
sty=str*sin(-xita)^2+stxi*cos(-xita)^2-2*tor*sin(-xita)*cos(-xita);
tau=-str*sin(-xita)*cos(-xita)+stxi*sin(-xita)*cos(-xita)+tor*(cos(-xita)^2-sin(-xita)^2);
end