R = 0.3;// define the radius of circle
L = 1.0;//

Point(1) = {L, -L, 0};//plot 7 points   first 4 points survey for square， last four for circle 
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};
Point(5) = {-L + R, -L, 0};
Point(6) = {-L, -L + R, 0};
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};

Circle(1) = {5, 4, 7};//draw the circle 
Circle(2) = {7, 4, 6};

Line(3) = {6, 3}; //connect the point to gennerate the mesh
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

Curve Loop(1) = {4, 7, 2, 3};// 用圆弧连接 定义平面
Plane Surface(1) = {1};

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;// define the length of line

Transfinite Surface{1};//参数化
Transfinite Surface{2};

Recombine Surface{1};//优化网格质量
Recombine Surface{2};

Mesh.ElementOrder = 1;//定义为线性网格
Mesh.Algorithm = 8;//设置网格算法为8

// EOF
