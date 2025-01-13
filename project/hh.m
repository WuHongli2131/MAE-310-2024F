
function [h_x, h_y] = hh (x, y, n_x, n_y)

L = 2;
r = @(x,y) sqrt((x+L/2)^2 + (y+L/2)^2);
th = @(x,y) atan2((y+L/2),(x+L/2));
Tx = 1e4;
R = 0.5;

xig_rr = @(x,y) Tx/2*(1-R^2/r(x,y)^2)+Tx/2*(1-4*R^2/r(x,y)^2+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_thth = @(x,y) Tx/2*(1+R^2/r(x,y)^2)-Tx/2*(1+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_rth = @(x,y) -Tx/2*(1+2*R^2/r(x,y)^2-3*R^4/r(x,y)^4)*sin(2*th(x,y));

xig_xx = @(x,y) xig_rr(x,y)*cos(-th(x,y))^2+xig_thth(x,y)*sin(-th(x,y))^2+2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_yy = @(x,y) xig_rr(x,y)*sin(-th(x,y))^2+xig_thth(x,y)*cos(-th(x,y))^2-2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_xy = @(x,y) -xig_rr(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_thth(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_rth(x,y)*(cos(-th(x,y))^2-sin(-th(x,y))^2);

xig_ij = zeros(2,2);
h = zeros(2,1);

vector = [n_x, n_y]';
xig_ij(1,1) = xig_xx(x,y);
xig_ij(1,2) = xig_xy(x,y);
xig_ij(2,1) = xig_xy(x,y);
xig_ij(2,2) = xig_yy(x,y);
for i = 1 : 2
    for j = 1 : 2
        h(i) = h(i) + xig_ij(i,j) * vector(j);
    end
end
h_x = h(1);
h_y = h(2);
end






function e_i = unit_vector(ii)
if ii == 1
    e_i = [1,0]';
elseif ii == 2
    e_i = [0,1]';
end
end


