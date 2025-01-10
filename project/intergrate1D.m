function val= intergrate1D(x1,y1,x2,y2,h1,h2)
%导入两个点的坐标，导入h 应该就可以计算出积分值
%线单元直接使用1阶积分  积分为 法向×形函数×h
[xi w]=Gauss(6,-1,1);
    x_l=zeros(6,1);
    y_l=zeros(6,1);
    d_dxi=zeros(6,1);
    dxi_d=x_l;
    for qua=1:6
        if x1==x2
            y_l(qua)= y1*PolyShape(1,1,xi(qua),0)+y2*PolyShape(1,2,xi(qua),0);
            d_dxi(qua)= y1*PolyShape(1,1,xi(qua),1)+y2*PolyShape(1,2,xi(qua),1);
        elseif y1==y2
            x_l(qua)=x1*PolyShape(1,1,xi(qua),0)+x2*PolyShape(1,2,xi(qua),0);
            d_dxi(qua)=x1*PolyShape(1,1,xi(qua),1)+x2*PolyShape(1,2,xi(qua),1);
        end
        dxi_d(qua)=1/d_dxi(qua);
    end
    f_el=0;
    for qua=1:6
        f_el=f_el+w(qua)* PolyShape(1, 1, xi(qua), 0)*h1*d_dxi(qua)+w(qua)* PolyShape(1, 2, xi(qua), 0)*h2*d_dxi(qua);
    end
    val=f_el;
end