function [xi, eta, weight] = Gauss2Dtri(n_int)
switch n_int
    case 6
        % preallocation
        xi = [0.659027622374092;0.659027622374092;0.23193368553031;0.23193368553031;
            0.109039009072877;0.109039009072877];
        eta=[0.23193368553031;0.109039009072877;0.659027622374092;0.109039009072877;
            0.659027622374092;0.23193368553031];
        % generate 1D rule
        weight=[1/6; 1/6; 1/6; 1/6; 1/6; 1/6];
    otherwise
        error("Illegal input!");
end
end
% EOF
%第三章
