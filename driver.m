%2D
clear all; clc;
%%整体思路：利用一个逻辑值来控制当前反向，然后在当前方向视为一个自由度二维问题，解出来后循环到下一个自由度，最后拼起来即可（当前版本不考虑3维）
%由于耦合，上述思路放弃了，转为直接使用书上思路
%%参数导入   需要用户修改的部分

E=1E9;%表示模量
mu=0.3;%泊松比
coe=E/(1-mu^2);%表系数
D = zeros(3);%D阵一步到胃
D(1,1)=coe;D(2,2)=D(1,1);D(1,2)=coe*-mu;D(2,1)=D(1,2);D(3,3)=coe*(1-mu)/2;
%%检验方程
% exact solution
dof=2;%自由度
dir=1;
v_x=0;%
v_y=0;
exact = @(x,y) (dir==1)*x*(1-x)*y*(1-y)+(dir==2)*x*(1-x)*y*(1-y);%dir表示当前方向，在前方加入循环即可将代码实现多次计算
exact_x = @(x,y) (dir==1)*(1-2*x)*y*(1-y)+(dir==2)*v_x;
exact_y = @(x,y) (dir==1)*x*(1-x)*(1-2*y)+(dir==2)*v_y;
f = @(x,y,d) (d==1)*((2*E*y*(y - 1))/(mu^2 - 1) - (E*(mu/2 - 1/2)...
    *((x - 1)*(y - 1) + x*y + 2*x*(x - 1) + x*(y - 1) + y*...
    (x - 1)))/(mu^2 - 1) + (E*mu*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(mu^2 - 1))...
    +(d==2)*((2*E*x*(x - 1))/(mu^2 - 1) - (E*(mu/2 - 1/2)*((x - 1)*(y - 1) + x*y + ...
    x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(mu^2 - 1) + ...
    (E*mu*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(mu^2 - 1)); % source term

%%四边形quadrature rule部分
% quadrature rule
n_int_xi  = 10;
n_int_eta = 10;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);%尝试修改 GAUSS 切换后生成高斯网格

%%生成网格部分，似乎已经被gmsh替代
% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = 6;               % number of elements in x-dir 划分单元格
n_el_y = 6;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements  总单元数  三角形中*2

n_np_x = n_el_x + 1;      % number of nodal points in x-dir  节点数
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points 总节点数

x_coor = zeros(n_np, 1);  %坐标 两者相同   划分为三角形时 每个矩形被划分为两个三角形 于是 n_el翻倍
y_coor = x_coor;            %但是n_np不变

hx = 1.0 / n_el_x;        % mesh size in x-dir  在三角形中hx需要*2  网格不均匀？
hy = 1.0 / n_el_y;        % mesh size in y-dir


% generate the nodal coordinates
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        index = (ny-1)*n_np_x + nx; % nodal index 似乎是用于鉴定节点是否正确的变量，无实际意义
        x_coor(index) = (nx-1) * hx;% 坐标  不用动
        y_coor(index) = (ny-1) * hy;
    end
end

%%IEN部分 这里可以加一个判断 跳过重新生成Ien
% IEN array

IEN = zeros(n_el, n_en); %element和其所连接的点的关系
for ex = 1 : n_el_x
    for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex; % element index  第几个元素
        IEN(ee, 1) = (ey-1) * n_np_x + ex;  %需要修改 原代码一次意外着一个方格，也就是一次可以描述两个三角形
        IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;%第一个三角形对应2*ee-1 第二个对应2*ee 再将节点补齐即可完成改写
        IEN(ee, 3) =  ey    * n_np_x + ex + 1;
        IEN(ee, 4) =  ey    * n_np_x + ex;
    end
end

% ID array
ID = zeros(n_np,2);
counter = 0;
for ny = 2 : n_np_y - 1%从2开始，因为边界没有方程
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;%为每一个点编号
        counter = counter + 1;
        ID(index,1) = counter;
        counter=counter+1;
        ID(index,2)=counter;
    end
end

n_eq = counter;%计算内部网格

%%k阵和f阵建立
% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq );%建立空阵 9还是很抽象  不用改
F = zeros(n_eq, 1);

%在这里理清一下思路 使用BDB矩阵表示K
%需要B阵，并且B阵与a和b无关
%需要向量e，不与a和b有关 和方向有关 脑子晕了这里   与i和j有关
%于是在一个单元内，定义出BDB 然后再与eiej循环？ 但是和Kele的关系？  Dij和eij绑定
%需要D阵，已定义
%额，还是在这个循环里面把误差一起弄了吧
%误差中k=1 但是

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );

    k_ele = zeros(2*n_en, 2*n_en); % element stiffness matrix 最后应该是一个三对角
    f_ele = zeros(2*n_en, 1);    % elemenlt oad vector

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;%xi表示可惜 eta表示伊塔 用于处理链式法则
        dy_dxi = 0.0; dy_deta = 0.0;
        Ba=zeros(3,2);%blackboard是一个3*2的矩阵 也就是B阵
        Bb=Ba;

        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll)); %取单个节点进行拟合
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));%一阶导拟合
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;%单节点的导数一起处理了
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;%雅可比行列式
        %因为可能要改三维。。这里ij写成循环形式
        for i=1:dof
            ei=(i==1)*[1,0]+(i==2)*[0,1];%还挺好用的这方式
            for aa = 1 : n_en%我又得去看书了，忘记原始公式了
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                %一阶导拟合 B阵在这里完成
                Ba(1,1)=Na_x;
                Ba(2,2)=Na_y;
                Ba(3,1)=Ba(2,2);
                Ba(3,2)=Ba(1,1);%单元内B阵完成
                pp=dof*(aa-1)+i;
                %出错点1？
                f_ele(pp) = f_ele(pp) + weight(ll) * detJ * f(x_l, y_l,i) * Na;%需要修改 似乎不用，i和j表示方向后就是看dir 但是两个方向是分开的 只需要最后存储的时候注意就可以了
                for j=1:dof
                    ej=(j==1)*[1,0]+(j==2)*[0,1];
                    for bb = 1 : n_en
                        qq=dof*(bb-1)+j;
                        Nb = Quad(bb, xi(ll), eta(ll));
                        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                        Bb(1,1)=Nb_x;
                        Bb(2,2)=Nb_y;
                        Bb(3,1)=Bb(2,2);
                        Bb(3,2)=Bb(1,1);
                        k_ele(pp,qq) = k_ele(pp,qq) + weight(ll) * detJ *ei*Ba'*D*Bb*ej';%确认正确
                    end % end of bb loop
                end % end of aa loop
            end % end of quadrature loop
        end
        %下面一段是原来ele和F的对应式子，需要补充边界条件就在这基础上修改
        for i=1:dof
            for aa = 1 : n_en
                PP = ID(IEN(ee,aa),i);
                if PP > 0  %比对1
                    F(PP) = F(PP) + f_ele(dof*(aa-1)+i);%算了不管了，反正我得到了K阵和F阵
                    for j=1:dof
                        for bb = 1 : n_en
                            QQ = ID(IEN(ee,bb),j);
                            if QQ > 0%比对2
                                K(PP, QQ) = K(PP, QQ) + k_ele(dof*(aa-1)+i, dof*(bb-1)+j);
                            else
                                % modify F with the boundary data
                                % here we do nothing because the boundary data g is zero or
                                % homogeneous
                            end
                        end
                    end
                end
            end
        end
        % F=F+(h,Na)+a(a,B)*g;
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
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");%为了proj做准备


% EOF
el2=0;
eh1=0;
for ee = 1 : n_el   %单元内划分
    x_ele = x_coor( IEN(ee, 1:n_en) );%节点内编号
    y_ele = y_coor( IEN(ee, 1:n_en) );%这里需要积分所在的精确的节点值
    %有disp，需要将disp和uh的单元对应起来   1：128 2：298 3：239 IEN阵
    uhp=disp(IEN(ee,:));%系数和节点对应起来  步骤1
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

            uh=uh+uhp(aa)*Na; %步骤2
            uh_x=uh_x+uhp(aa)*Na_x;%存疑项
            uh_y=uh_y+uhp(aa)*Na_y;

            %计算uh
        end % end of aa loop
        %到这里暂存了一个节点的uh uh_x和uh_y   并且在quadrature rule 内部 这里完成积分
        e=uh-exact(x_l,y_l);
        e_x=uh_x-exact_x(x_l,y_l);%可能需要两个自由度分别表示 但是耦合项？
        e_y=uh_y-exact_y(x_l,y_l);
        el2=el2+weight(ll)*detJ*(e^2);
        eh1=eh1+weight(ll)*detJ*(e_x^2+e_y^2);
    end % end of quadrature loop
end
yy(i)=sqrt(el2);
yd(i)=sqrt(eh1);
xx(i)=hx;
