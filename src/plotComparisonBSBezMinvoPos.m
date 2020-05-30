close all; clc; clear;


set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


syms t real

addpath(genpath('./utils'));

interv=[0,1];


T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';



%% Plot volume MINVO
A=getSolutionA(3,"01");

figure;
subplot(2,2,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 2000])



v0=[1.1    -1/sqrt(3)   0]';
v1=[-0.5   -1/sqrt(3)   0.3]';
v2=[0     2/sqrt(3)   0.8]';
v3=[0.3     0           4/sqrt(6)]';

V=[v0 v1 v2 v3];
vx=V(1,:)';
vy=V(2,:)';
vz=V(3,:)';


pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

P=[pol_x'; pol_y'; pol_z'];

volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.05);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);




%% Plot volume Bezier
subplot(2,2,2);hold on
A=computeMatrixForBezier(3,"01");
% pol_x=A'*vx;
% pol_y=A'*vy;
% pol_z=A'*vz;
volumen_bezier=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.05);

fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);

 
%% Plot volume B-Spline
subplot(2,2,[3,4]);hold on

Mbs=(1/6)*[1 -3 3 -1
                4  0  -6  3
                1  3  3  -3
               0  0  0  1];

A=[Mbs(:,4), Mbs(:,3), Mbs(:,2), Mbs(:,1)]; 


V=P*inv(A);   %P=VA

p = 3;
n = 4; %4 for degree 3
m=n+p+1;
num_of_intervals=m-2*p;
knots = [0     0     0     0  0.143 0.286 0.429 0.571 0.714 0.857     1     1     1     1];  % knot vector
       
       
cPoints=[  0  0.143  0.429;
          0   0.143  0.429;
          1   1.14   1.43];

 

 
cPoints=[cPoints V];


tmp=[    9.07   9.07   9.07;
       -1.07  -1.07  -1.07;
        -1.2   -1.2   -1.2];
 
cPoints=[cPoints tmp];
    
increm=0.0001;

num_points=200;
points_bs = bspline_deboor(n,knots,cPoints,num_points );%0:increm:1


points_bs_cropped=points_bs(:,0.4*num_points:0.6*num_points);


                                


volumen_bs=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.3);
subplot(2,2,[3,4]);hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,[min(knots),max(knots)],'r','LineWidth',3);

plot3(points_bs_cropped(1,:),points_bs_cropped(2,:),points_bs_cropped(3,:),'-b')

subplot(2,2,1);hold on
plot3(points_bs_cropped(1,:),points_bs_cropped(2,:),points_bs_cropped(3,:),'-b')

subplot(2,2,2);hold on
plot3(points_bs_cropped(1,:),points_bs_cropped(2,:),points_bs_cropped(3,:),'-b')


view1=80;
view2=30;

subplot(2,2,1);
title(['\textbf{MINVO, Volume=',num2str(volumen_minvo,4),'  $u^3$}'])
view(view1, view2)
axis equal;
xlim([-1,2]);
ylim([-2,2]);
arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);
axis off

subplot(2,2,2);
title(['\textbf{B\''ezier, Volume=',num2str(volumen_bezier,4),'  $u^3$}'])
view(view1, view2)
axis equal;
arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

axis off
% ylim([-1.5,1.5]);

subplot(2,2,[3,4]);
title(['\textbf{B-Spline, Volume=',num2str(volumen_bs,4),'  $u^3$}'])
view(view1, view2)
axis equal;
arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
axis off
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);
% ylim([-1.5,1.5]);

% exportAsSvg(gcf,'imgs/comparisonBsBezierMinvo_matlab')