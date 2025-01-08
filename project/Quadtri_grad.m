function [val_xi, val_eta] = Quadtri_grad(aa, xi, eta)

if aa == 1  %导数居然与位置无关
    val_xi  = -1;
    val_eta = -1;
elseif aa == 2
    val_xi  =  1;
    val_eta = 0;
elseif aa == 3
    val_xi  = 0;
    val_eta = 1;
else
    error('Error: value of a should be 1,2,3, or 4.');
end

% EOF