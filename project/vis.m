close all;
[X, Y] = meshgrid(0 : hx : 1, 0 : hy : 1);
hold on;
Z = reshape(disp(:,2), n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp;
axis equal;