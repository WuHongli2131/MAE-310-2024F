clear all; clc;  % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh

pp   = 2;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes 3
err=zeros(8,1);        % 创建一个数列储存误差
eh1=err;               % first order error
for n_el=2:2:16      % number of elements 使得函数hh从2到16循环
    n_np = n_el * pp + 1;  % number of nodal points  5
    n_eq = n_np - 1;       % number of equations  4
    n_int = 10;

    hh = 1.0 / (n_np - 1); % space between two adjacent nodes
    x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

    IEN = zeros(n_el, n_en);

    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end

    % Setup the ID array for the problem
    ID = 1 : n_np;
    ID(end) = 0;

    % Setup the quadrature rule
    [xi, weight] = Gauss(n_int, -1, 1);

    % allocate the stiffness matrix
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);% 快速三角阵
    F = zeros(n_eq, 1);

    % Assembly of the stiffness matrix and load vector  一个单元的
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
        f_ele = zeros(n_en, 1);    % allocate a zero element load vector

        x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

        % quadrature loop
        for qua = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                end
            end
        end

        % Assembly of the matrix and vector based on the ID or LM data
        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                    end
                end
            end
        end
    end

    % ee = 1 F = NA(0)xh
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

    % Solve Kd = F equation
    d_temp = K \ F;

    disp = [d_temp; g];

    % Postprocessing: visualization
    %plot(x_coor, disp, '--r','LineWidth',3);
    %x_sam = 0 : 0.01 : 1;
    %y_sam = x_sam.^5;
    %hold on;
    %plot(x_sam, y_sam, '-k', 'LineWidth', 3);

    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;

    x_sam = zeros(n_el * n_sam + 1, 1);% coordinate
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sa = x_sam; % store the numerical solution value at sampling pts
    ud_sam= x_sam;
    yd_sam= x_sam;
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        if ee == n_el
            n_sam_end = n_sam+1;
        else
            n_sam_end = n_sam;
        end

        for ll = 1 : n_sam_end
            x_l = 0.0;
            u_l = 0.0;
            u_sd=0;    %d means dot
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_sd = u_sd + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1);
            end
            x_sam( (ee-1)*n_sam + ll ) = x_l;
            u_sam( (ee-1)*n_sam + ll ) = u_l;
            ud_sam(( (ee-1)*n_sam + ll ))=u_sd;
            y_sam( (ee-1)*n_sam + ll ) = x_l^5;
            yd_sam( (ee-1)*n_sam + ll ) =5*x_l^4;
        end
    end
    % calculate the integration
    u_f=zeros(n_el,1);
    u_fd=u_f;%means final
    y_f=u_f;
    yd_f=u_f;
    for mm=1:n_el
        [ui, weight2] = Gauss(n_int, x_coor((mm-1)*pp+1), x_coor(mm*pp+1));%更换积分区间 一阶表示ui
        x1sam=ui;
        for ii=1:n_int  %小区间积分
            u_s=0;
            u_sd=0;
            y_s=0;
            y_sd=0;
            for aa = 1 : n_en
                u_s = u_s + ui(aa) * PolyShape(pp, aa, x1sam(ii), 0);%积分所需u
                u_sd = u_sd + ui(aa) * PolyShape(pp, aa, x1sam(ii), 1);
                y_s = y_s + ui(aa) * PolyShape(pp, aa, y_sam(ii), 0);%积分所需u
                y_sd = y_sd + ui(aa) * PolyShape(pp, aa, yd_sam(ii), 1);
            end
            u_f(mm)=u_s;
            u_fd(mm)=u_sd;
            y_f(mm)=y_s;
            y_fd(mm)=y_sd;
        end             %得到x u u一阶导
        %quardratic将积分弄出，最好是存到某一数组里面
        %求和，利用quadratic把y弄出了，最后得到e
        %待实现
    end
    [yi,weight3]=Gauss(n_el, 0, 1);
    integ1=0;          %临时积分变量
    integ2=0; %什么天才会把临时变量放到最外面啊？？？？？！！！@@#￥%……
    integ_1=0;
    integ_2=0;
    for nn=1:length(yi)
        integ1=integ1+weight3(nn)*(u_f(nn)-y_f(nn)).^2;  %发现bug，nn含义不同，导致积分点选错  需要用nn表示ui（找不到明显规律）
        integ2=integ2+weight3(nn)*u_f(nn).^2;
        integ_1=integ_1+weight3(nn)*(u_fd(nn)-y_fd(nn)).^2;
        integ_2=integ_2+weight3(nn)*(u_fd(nn)).^2;
    end

    err(n_el/2)=integ1/integ2;
    eh1(n_el/2)=integ_1/integ_2;
end
%





























% EOF