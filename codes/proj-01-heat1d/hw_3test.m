clear; clc;  % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

exact = @(x) x.^5;
exact_x = @(x) 5 * x.^4;
% Setup the mesh
err=zeros(8,1);        % 创建一个数列储存误差  %第二次天才
eh1=err;               % first order error

for pp   = 1:3         % polynomial degree
    n_en = pp + 1;         % number of element or local nodes 3
    xx=zeros(8,1);          % ready for plot  in the loop
    yy=xx;
    yd=xx;
    for n_el=2:2:16   % number of elements 使得函数hh从2到16循环
        n_np = n_el * pp + 1;  % number of nodal points  5
        n_eq = n_np - 1;       % number of equations  4
        n_int = 20;

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
            % 关于映射  首先映射到一个qua坐标系下，qua表示编号，而x_l表示的是新系，1 在新系下换元可以求导用于积分。
            for qua = 1 : n_int
                dx_dxi = 0.0;
                x_l = 0.0;
                for aa = 1 : n_en
                    x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);% 可惜
                    dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);%d可惜dx
                end
                dxi_dx = 1.0 / dx_dxi;%dxd可惜

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
        % 
        % n_sam = 20;
        % xi_sam = -1 : (2/n_sam) : 1;
        % 
        % x_qua = zeros(n_el * n_sam + 1, 1);% coordinate
        % y_sam = x_qua; % store the exact solution value at sampling points
        % u_sam = x_qua; % store the numerical solution value at sampling pts
        % ud_sam= x_qua;
        % yd_sam= x_qua;
        % for ee = 1 : n_el
        %     x_ele = x_coor( IEN(ee, :) );
        %     u_ele = disp( IEN(ee, :) );
        % 
        %     if ee == n_el
        %         n_sam_end = n_sam+1;
        %     else
        %         n_sam_end = n_sam;
        %     end
        % 
        %     for ll = 1 : n_sam_end
        %         x_l = 0.0;
        %         u_l = 0.0;
        %         u_sd=0;    %d means dot
        %         for aa = 1 : n_en
        %             x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);%精确值*形函数
        %             u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
        %             u_sd = u_sd + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1);
        %         end
        %         x_qua( (ee-1)*n_sam + ll ) = x_l;%coor
        %         u_sam( (ee-1)*n_sam + ll ) = u_l;%uh
        %         ud_sam(( (ee-1)*n_sam + ll ))=u_sd;%uh、
        %         y_sam( (ee-1)*n_sam + ll ) = x_l^5;%y
        %         yd_sam( (ee-1)*n_sam + ll ) =5*x_l^4;%yh、
        %     end
        % end
        % calculate the integration
        u_f=zeros(n_el,1);
        u_fd=u_f;%means final
        y_f=u_f;
        yd_f=u_f;

        integ1=0;          %临时积分变量
        integ2=0;           %什么天才会把临时变量放到最外面啊？？？？？！！！@@#￥%……
        integ_1=0;
        integ_2=0;
        %伪代码：首先已知节点精确，可用polyshape表示内部的点，得到高斯积分所需值，然后加权求和
        xh=zeros(n_int,1);
        uh=xh;
        uhd=xh;
        %积重难反，嵌套太多，重写了。
        %第一步：需要高斯积分点的uh  选择复制之前的代码 将其中的n_sam替换成我需要的ni
     L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;


        for ee = 1 : n_el
            x_ele = x_coor( IEN(ee, :) );
            u_ele = disp( IEN(ee, :) );
            for ll = 1 : n_int
                x_l = 0.0;
                dx_dxi=0;
                u_l = 0.0;
                u_sd=0;    %d means dot
                for aa=1:n_en%       bug
                    x_l=x_l+x_ele(aa)*PolyShape(pp,aa,xi(ll),0);%这里就得到需要的u和uh了
                    dx_dxi=dx_dxi+x_ele(aa)*PolyShape(pp,aa,xi(ll),1);
                    u_l=u_l+u_ele(aa)*PolyShape(pp,aa,xi(ll),0);
                    u_sd=u_sd+u_ele(aa)*PolyShape(pp,aa,xi(ll),1);%对kexi求导 bug
                end
                dxi_dx=1/dx_dxi;
               %高斯积分
                integ1= integ1+weight(ll)*(u_l-exact(x_l))^2.*dx_dxi;%ele内部积分得到完全积分
                integ2= integ2+weight(ll)*exact(x_l)^2.*dx_dxi;
                integ_1=integ_1+weight(ll)* ( u_sd * dxi_dx - exact_x(x_l) )^2*dx_dxi;
                integ_2=integ_2+weight(ll)*exact_x(x_l)^2*dx_dxi;
                %
                % x_qua( (ee-1)*n_int + ll ) = x_l;%coor
                % u_q( (ee-1)*n_int + ll ) = u_l;%uh
                % ud_q(( (ee-1)*n_int + ll ))=u_sd;%uh、
                % y_q( (ee-1)*n_int + ll ) = x_l^5;%y
                % yd_q( (ee-1)*n_int + ll ) =5*x_l^4;%yh、
            end
        upper=sqrt(integ1);
        lower=sqrt(integ2);
        upperd=sqrt(integ_1);
        lowerd=sqrt(integ_2);
        end
       
        err(n_el/2)=upper/lower;
        eh1(n_el/2)=upperd/lowerd;
        xx(n_el/2)=log(1/n_el);
        yy(n_el/2)=log(err(n_el/2));
        yd(n_el/2)=log(eh1(n_el/2));

        % for mm=1:n_int %问题都在这段积分里
        %     %u在结点处不精确
        %
        %     for qua = 1 : n_el         %学习quardratic结果
        %         x_ele = x_coor( IEN(qua, :) );
        %         u_ele = disp( IEN(qua, :) );
        %         if qua == n_el
        %         n_sam_end = n_sam+1;
        %     else
        %         n_sam_end = n_sam;
        %     end
        %         dx_dxi = 0.0;
        %         x_l = 0.0;
        %         u_s=0;
        %         u_sd=0;
        %
        %         for aa = 1 : n_en%得到quadratic点处所需的uh值   %这里积分异常  节点不精确，因为用的是u——sam进行了2此近似。
        %             x_l    = x_l + x_ele(aa) * PolyShape(pp, aa, xi(mm), 0);
        %             u_s   = u_s+ u_ele(aa) * PolyShape(pp, aa, xi(mm), 0);
        %             dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(mm), 1);%这是dxd可惜ok
        %             u_sd  = u_sd  + u_ele(aa) * PolyShape(pp, aa, xi(mm), 1);
        %         end
        %         dxi_dx = 1.0 / dx_dxi;
        %         xh(qua)=x_l;%记录quardradic  这里已经是节点精确值了
        %         uh(qua)=u_s;
        %         uhd(qua)=u_sd;
        %         % u_f(mm) = u_f(mm) + weight(qua) * (uh(qua)-xh(qua).^5).^2 * dx_dxi;%容易出错的循环嵌套
        %         % u_fd(mm)=u_fd(mm)+weight(qua)*(uhd(qua)-5*xh(qua).^4).^2*dx_dxi;%相当于积分 次数是nel次加权求和
        %         % y_f(mm)=y_f(mm)+weight(qua)*uh(qua).^2*dx_dxi;
        %         % yd_f(mm)=yd_f(mm)+weight(qua)*uhd(qua).^2*dx_dxi;
        %     end%循环出来后u_f（mm）操作有问题
        %
        %         integ1=integ1+weight(mm) * (uh(mm)-xh(mm).^5).^2*dxi_dx;  %发现bug，nn含义不同，导致积分点选错  需要用nn表示ui（找不到明显规律）
        %         integ2=integ2+weight(mm)*(uhd(mm)-5*xh(mm).^4).^2*dxi_dx;
        %         integ_1=integ_1+weight(mm)*(xh(mm).^5).^2*dxi_dx;%使用y_sam和yd_sam直接得到积分
        %         integ_2= integ_2+weight(mm)*(5*xh(mm)).^4.^2*dxi_dx;
        % end
        %     err(n_el/2,pp-1)=sqrt(integ_1)/sqrt(integ1);
        %     eh1(n_el/2,pp-1)=sqrt(integ_2)/sqrt(integ2);%
        %     xx(n_el/2)=log(hh);
        %     yy(n_el/2)=log(err(n_el/2,pp-1));
        %     yd(n_el/2)=log(eh1(n_el/2,pp-1));
        % end
        % figure;
        % plot(xx,yy, '-r','LineWidth',3);
        % hold on;
        % plot(xx,yd, '-k','LineWidth',3);
    end
   
    figure;
    plot(xx,yy, '-r','LineWidth',3);
    hold on;
    plot(xx,yd, '-k','LineWidth',3);%画出图像仍然有问题
end
% EOF