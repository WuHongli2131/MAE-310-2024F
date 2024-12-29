function val = Quadtri(aa, xi, eta)

if aa == 1
    val = 1-xi-eta;
elseif aa == 2
    val = xi;
elseif aa == 3
    val = eta;
else
    error('Error: value of a should be 1,2,3, or 4.');
end

% EOF