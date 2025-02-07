%2D
clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int= 6;
[xi, eta, weight] = Gauss2Dtri(n_int);%尝试修改 GAUSS 切换后生成高斯网格

% mesh generation
n_en= 3;% number of nodes in an element  %3个节点
xx=zeros(7,1);
yy=xx;
yd=xx;
i=0;
for n_el_x = 2:2:14    % number of elements in x-dir 划分单元格 误差修改这里以修改h
    i=i+1;
    n_el_y = n_el_x;               % number of elements in y-dir
    n_el   =2* n_el_x * n_el_y; % total number of elements  总单元数  三角形中*2

    n_np_x = n_el_x + 1;      % number of nodal points in x-dir  节点数
    n_np_y = n_el_y + 1;      % number of nodal points in y-dir
    n_np   = n_np_x * n_np_y; % total number of nodal points 总节点数

    x_coor = zeros(n_np, 1);  %坐标 两者相同   划分为三角形时 每个矩形被划分为两个三角形 于是 n_el翻倍
    y_coor = x_coor;            %但是n_np不变

    hx = 1.0 / n_el_x;        % mesh size in x-dir  在三角形中hx需要*2 这里不用×2，因为网格只是被分成两部分，没有被改变
    hy = 1.0 / n_el_y;        % mesh size in y-dir

    % generate the nodal coordinates
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx; % nodal index 似乎是用于鉴定节点是否正确的变量，无实际意义
            x_coor(index) = (nx-1) * hx;% 坐标  不用动
            y_coor(index) = (ny-1) * hy;
        end
    end

    % IEN array
    IEN = zeros(n_el, n_en); %element和其所连接的点的关系
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index  第几个元素
            IEN(2*ee-1, 1) = (ey-1) * n_np_x + ex;  %需要修改 原代码一次意外着一个方格，也就是一次可以描述两个三角形
            IEN(2*ee-1, 2) = (ey-1) * n_np_x + ex + 1;%第一个三角形对应2*ee-1 第二个对应2*ee 再将节点补齐即可完成改写

            IEN(2*ee-1, 3) =  ey    * n_np_x + ex;%第一个三角形 124
            %第二个三角形 234
            IEN(2*ee, 1) = (ey-1) * n_np_x + ex + 1;
            IEN(2*ee, 2) =  ey    * n_np_x + ex + 1;
            IEN(2*ee, 3) =  ey    * n_np_x + ex;
        end%已验证没有问题
    end

    % ID array
    ID = zeros(n_np,1);
    counter = 0;
    for ny = 2 : n_np_y - 1%从2开始，因为边界没有方程
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;%为每一个点编号 确认其有方程 方程数似乎只和节点数有关 不需要修改？
            counter = counter + 1;%验证完毕
            ID(index) = counter;
        end
    end

    n_eq = counter;%计算内部网格

    LM = ID(IEN);%略显抽象，只知道结果是把方程映射到了节点上  验证完毕

    % allocate the stiffness matrix and load vector
    K = spalloc(n_eq, n_eq, 9 * n_eq);%建立空阵 9还是很抽象  不用改
    F = zeros(n_eq, 1);

    % loop over element to assembly the matrix and vector
    for ee = 1 : n_el   %单元内划分
        x_ele = x_coor( IEN(ee, 1:n_en) );%节点内编号
        y_ele = y_coor( IEN(ee, 1:n_en) );

        k_ele = zeros(n_en, n_en); % element stiffness matrix
        f_ele = zeros(n_en, 1);    % element load vector

        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;%xi表示可惜 eta表示伊塔 用于处理链式法则
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quadtri(aa, xi(ll), eta(ll)); %取单个节点进行拟合
                y_l = y_l + y_ele(aa) * Quadtri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quadtri_grad(aa, xi(ll), eta(ll));%一阶导拟合
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;%单节点的导数一起处理了
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;%雅可比行列式  到这里都没问题

            for aa = 1 : n_en%我又得去看书了，忘记原始公式了 拟合的准则似乎没变 变得只有形函数
                Na = Quadtri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quadtri_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

                for bb = 1 : n_en
                    Nb = Quadtri(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quadtri_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end % end of bb loop
            end % end of aa loop
        end % end of quadrature loop

        for aa = 1 : n_en
            PP = LM(ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(aa);%算了不管了，反正我得到了K阵和F阵

                for bb = 1 : n_en
                    QQ = LM(ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % modify F with the boundary data
                        % here we do nothing because the boundary data g is zero or
                        % homogeneous
                    end
                end
            end
        end
    end

    % solve the stiffness matrix
    dn = K \ F;

    % insert dn back into the vector for all nodes
    disp = zeros(n_np, 1); %得到拟合项系数

    for ii = 1 : n_np
        index = ID(ii);
        if index > 0
            disp(ii) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end% save the solution vector and number of elements to disp with name
    % HEAT.mat
    save("HEAT", "disp", "n_el_x", "n_el_y");%为了proj做准备

    % EOF
    %第一遍结果貌似有bug 出来和原来的图像不一致

    %%下面是关于err的代码
    %从书上了解到err应该满足一个关系，在三角形中形函数的k=1，m=1或0 k+1-m=2或1
    %m=0时应该是h^2关系， m=1时应该是h的关系取log看斜率即可
    %需要积分uh-u的平方 总之前面积分部分没有问题，新开一个积分式子就行了
    %具体动手写，你也是个天才，神特么用键盘推积分公式
    %最后推出来发现行列式在笔记本上 服了
    el2=0;
    eh1=0;
    for ee = 1 : n_el   %单元内划分
        x_ele = x_coor( IEN(ee, 1:n_en) );%节点内编号
        y_ele = y_coor( IEN(ee, 1:n_en) );%这里需要积分所在的精确的节点值
        %有disp，需要将disp和uh的单元对应起来   1：128 2：298 3：239 IEN阵
        uhp=disp(IEN(ee,:));%系数和节点对应起来
        for ll = 1 : n_int  %开始积分 重复上述，大部分不用改
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;%xi表示可惜 eta表示伊塔 用于处理链式法则
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quadtri(aa, xi(ll), eta(ll)); %取单个节点进行拟合
                y_l = y_l + y_ele(aa) * Quadtri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quadtri_grad(aa, xi(ll), eta(ll));%一阶导拟合
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;%单节点的导数一起处理了
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;%雅可比行列式  到这里都没问题

            uh=0; uh_x=0;uh_y=0;u=0;
            for aa = 1 : n_en%我又得去看书了，忘记原始公式了 拟合的准则似乎没变 变得只有形函数
                Na = Quadtri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quadtri_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;%笔记公式
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                % 这里不用b，只用计算一次积分
                % f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
                %
                % for bb = 1 : n_en
                %   Nb = Quad(bb, xi(ll), eta(ll));
                %   [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                %   Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                %   Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                %
                %   k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                %
                % end % end of bb loop
                uh=uh+uhp(aa)*Na;
                uh_x=uh_x+uhp(aa)*Na_x;
                uh_y=uh_y+uhp(aa)*Na_y;

                %计算uh
            end % end of aa loop
            %到这里暂存了一个节点的uh uh_x和uh_y   并且在quadrature rule 内部 这里完成积分
            e=uh-exact(x_l,y_l);
            e_x=uh_x-exact_x(x_l,y_l);
            e_y=uh_y-exact_y(x_l,y_l);
            el2=el2+weight(ll)*detJ*(e^2);
            eh1=eh1+weight(ll)*detJ*(e_x^2+e_y^2);
        end % end of quadrature loop
    end
    yy(i)=sqrt(el2);
    yd(i)=sqrt(eh1);
    xx(i)=hx;
end
figure;
plot(log(xx), log(yy), '-r','LineWidth',3);%出来函数图像很奇怪 不知道哪里出问题 拟合结果为4和3
hold on;
plot(log(xx),log(yd),'b');