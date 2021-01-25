%T is the 4x4 pose matrix: [R t; 0 1]
function plotAxesArrowsTnoColor(length, w_T_b)

w_t_b=w_T_b(1:3,4);

%Let D be a point along x_body, E along y_body,...

b_D=[length 0 0 1]';
b_E=[0 length 0 1]';
b_F=[0 0 length 1]';

w_D=w_T_b*b_D;
w_E=w_T_b*b_E;
w_F=w_T_b*b_F;



arrow3d(w_t_b',w_D(1:3)',20,'cylinder',[0.2,0.1]);
arrow3d(w_t_b',w_E(1:3)',20,'cylinder',[0.2,0.1]);
arrow3d(w_t_b',w_F(1:3)',20,'cylinder',[0.2,0.1]);
end