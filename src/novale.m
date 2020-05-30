
close all; clear; clc;

addpath(genpath('./utils'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


name_figure='imgs/plots_basis';

syms t real
interv=[-1,1];
figure;
n_rows=6;
n_cols=4;


T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';

%%GEOMETRIC INTERPRETATION OF LAMBDA_I
figure;hold on
% subplot(1,2,1);
set(gcf, 'Position',  [500, 500, 3000, 1000])

A=getSolutionA(3,"m11");
view1=150;
view2=30;

v0=[0 0 0]';
v1=[1 0 0]';
v2=[0 1 0]';
v3=[0 0 1]';

V=[v0 v1 v2 v3];
vx=V(1,:)';
vy=V(2,:)';
vz=V(3,:)';

pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

pol_x=[0 0 1 0]';
pol_y=[0 1 0 0]';
pol_z=[1 0 0 0]';

P=[pol_x'; pol_y'; pol_z']

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
 axis equal;
% ylim([-1.5,1.5]);

% arrow3d([0 0 0],[0 0 1.5],10,'cylinder',[0.2,0.1]);
% arrow3d([0 0 0],[0 1.5 0],10,'cylinder',[0.2,0.1]);
% arrow3d([0 0 0],[1.5 0 0],10,'cylinder',[0.2,0.1]);


AinvT=inv(A)'
% figure;
x=AinvT(:,1)
y=AinvT(:,2)
z=AinvT(:,3)
[k1,av1] = convhull(x,y,z);
trisurf(k1,x,y,z,'FaceColor','cyan')
axis equal

% plotSphere(position, radius, color)


% Plot the parallelepiped
% pt0=[pol_x(end) pol_y(end) pol_z(end)]';
% p0=[0 0 0]';
% p1=P(:,1);%[]pol_x(1:end-1);
% p2=P(:,2);%pol_y(1:end-1);
% p3=P(:,3);%pol_z(1:end-1);
% p4=p1+p2;
% p5=p2+p3;
% p6=p1+p3;
% p7=p1+p2+p3;
% pointsParallep=[p0+pt0 p1+pt0 p2+pt0 p3+pt0 p4+pt0 p5+pt0 p6+pt0 p7+pt0];
% 
% [k1,av1] = convhull(pointsParallep(1,:),pointsParallep(2,:),pointsParallep(3,:));
% 
% trisurf(k1,pointsParallep(1,:),pointsParallep(2,:),pointsParallep(3,:),'FaceColor','cyan')


%% And now on that figure, plot the solution on the base



%% RESULT for 2D for a given simplex

v0=v1(1:2)
v1=v2(1:2)
v2=v3(1:2)

A=getSolutionA(2,"m11");


vx=[v1(1)  v2(1)  v3(1)]';
vy=[v1(2)  v2(2)  v3(2)]';

V=[v0 v1 v2];


pol_x=pol_x(2:end);
pol_y=pol_y(2:end);
P=[pol_x'; pol_y'];

V=P*inv(A);

[k,av] = convhull(V');
fill3(V(1,k),V(2,k),zeros(size(V(1,k))),'b','LineWidth',1)
alpha 0.2
axis equal
%%
% ylim([0 2.0])
% pol_x=A'*vx;
% pol_y=A'*vy;
% 
% plot(0,0,'*')
fplot(pol_x'*T2,pol_y'*T2,interv,'b','LineWidth',3);
% ylim([-0.1,inf])

for(i=1:size(A,1))
    roots_poly=real(roots(A(i,:))); %real to avoid numerical approximations
    x=subs(pol_x'*T2,t,roots_poly(1));
    y=subs(pol_y'*T2,t,roots_poly(1));
    plot3(x,y,0.0,'o')
end


%%
A1=rand(3,1);
A2=rand(3,3);
V=rand(2,3);

P=V*[A1 A2]

V*A1

V*A2


%%
% clc
% pol_x=[0 2 0 4]';
% pol_y=[0 1 0 8]';
% pol_z=[1 2 0 1]';
% 
% P=[pol_x'; pol_y'; pol_z'];
% V=P*inv(A)
% 
% det([V' ones(4,1)]) %Check if the points are coplanar


