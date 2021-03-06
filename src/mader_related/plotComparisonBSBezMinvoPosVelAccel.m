% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


close all; clc; clear;
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

set(0, 'DefaultFigureRenderer', 'painters'); %faster than opengl

syms t real

addpath(genpath('../solutions'));
addpath(genpath('../utils'));

interv=[0,1];

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';

v0=[-1/sqrt(3)    1.1   0]';
v1=[-0.5   -1/sqrt(3)   0.3]';
v2=[0     2/sqrt(3)   0.8]';
v3=[1.6     0           4/sqrt(6)]';

%For the MINVO paper use this:
% v0=[1.1    -1/sqrt(3)   0]';
% v1=[-0.5   -1/sqrt(3)   0.3]';
% v2=[0     2/sqrt(3)   0.8]';
% v3=[0.3     0           4/sqrt(6)]';

V=[v0 v1 v2 v3];
vx=V(1,:)'; vy=V(2,:)'; vz=V(3,:)';
A=getA_MV(3,interv);
P=V*A;
pol_x=P(1,:)'; pol_y=P(2,:)'; pol_z=P(3,:)';
v_x=polyder(pol_x)'; v_y=polyder(pol_y)'; v_z=polyder(pol_z)';
a_x=polyder(v_x)'; a_y=polyder(v_y)'; a_z=polyder(v_z)';


%% Plot volume MINVO
figure; hold on; set(gcf, 'Position',  [500, 500, 3000, 2000])
subplot(3,3,1);hold on;title('Position');xlabel('x'); ylabel('y'); zlabel('z') %set(gcf, 'Position',  [500, 500, 3000, 2000])

%%Position
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.07);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
xlabel('x'); ylabel('y'); zlabel('z')
title(['\textbf{Pos, vol=',num2str(volumen_minvo,4),'  $u_p^3$}'])
%%
%%Velocity
%figure; hold on;
subplot(3,3,2); hold on; title('Velocity'); xlabel('vx'); ylabel('vy'); zlabel('vz')
A=getA_MV(2,interv);
area_minvo= plot_plane_convex_hull(v_x, v_y, v_z, A, 'g', 0.5);
fplot3(v_x'*T2,v_y'*T2,v_z'*T2,interv,'r','LineWidth',3);
title(['\textbf{Vel, area=',num2str(area_minvo,4),'  $u_v^2$}'])
%%
%%Acceleration
subplot(3,3,3); hold on; title('Acceleration'); xlabel('ax'); ylabel('ay'); zlabel('az')
A=getA_MV(1,interv);
long_minvo=plot_line_convex_hull(a_x,a_y,a_z,A,'g',1.5);
fplot3(a_x'*T1,a_y'*T1,a_z'*T1,interv,'r','LineWidth',3);
title(['\textbf{Accel, long=',num2str(long_minvo,4),'  $u_a$}'])

%% Plot volume Bezier
%figure; hold on;
subplot(3,3,4);hold on;xlabel('x'); ylabel('y'); zlabel('z')
A=getA_Be(3,interv);
volumen_bezier=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.07);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
title(['\textbf{Pos, vol=',num2str(volumen_bezier,4),'  $u_p^3$}'])
%%
%%Velocity
subplot(3,3,5); hold on;  xlabel('vx'); ylabel('vy'); zlabel('vz')
A=getA_Be(2,interv);
area_bezier= plot_plane_convex_hull(v_x, v_y, v_z, A, 'b', 0.5);
fplot3(v_x'*T2,v_y'*T2,v_z'*T2,interv,'r','LineWidth',3);
title(['\textbf{Vel, area=',num2str(area_bezier,4),'  $u_v^2$}'])

%%Acceleration
subplot(3,3,6); hold on; xlabel('ax'); ylabel('ay'); zlabel('az')
A=getA_Be(1,interv);
long_bezier= plot_line_convex_hull(a_x,a_y,a_z,A,'b',1.5);
fplot3(a_x'*T1,a_y'*T1,a_z'*T1,interv,'r','LineWidth',3);
title(['\textbf{Accel, long=',num2str(long_bezier,4),'  $u_a$}'])
 
%% Plot volume B-Spline
subplot(3,3,7); hold on; title('Position'); xlabel('x'); ylabel('y'); zlabel('z')

A=getA_BS(3,interv);

V=P*inv(A);   %P=VA

p = 3;
deltaT=1/7;
knots = [0     0     0     0  deltaT*[1:16]   19*deltaT*ones(1,4)];  % knot vector
cPoints=[  0  0.143  0.429 0.5  1.1  7.3  1.6  1.9;
          0   0.143  0.429 0.7  1.3  7.5  1.7 1.5;
          1   1.14   1.43  0.9   1.5  7.2  1.8  1.6];
cPoints=[cPoints V];
tmp=[  0.3  0.8  0.4  0.6  0.8  9.07   9.07   9.07;
       -1.2 1.1 0.4  2.3  2.5 -1.07  -1.07  -1.07;
       -1.4  1.5  2.6  2.4  4.5 -1.2   -1.2   -1.2];
cPoints=[cPoints tmp];

                            
volumen_bs=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.5);

fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,[0,1],'r','LineWidth',3);
title(['\textbf{Pos, vol=',num2str(volumen_bs,4),'  $u_p^3$}'])

%%Velocity
subplot(3,3,8); hold on; title('Velocity'); xlabel('vx'); ylabel('vy'); zlabel('vz')
A=getA_BS(2,interv);
area_bs= plot_plane_convex_hull(v_x, v_y, v_z, A, 'y', 0.5);
fplot3(v_x'*T2,v_y'*T2,v_z'*T2,interv,'r','LineWidth',3);
title(['\textbf{Vel, area=',num2str(area_bs,4),'  $u_v^2$}'])
%% 
%% Acceleration
%figure; hold on;
subplot(3,3,9);  hold on; xlabel('ax'); ylabel('ay'); zlabel('az')
A=getA_BS(1,interv);
long_bs= plot_line_convex_hull(a_x,a_y,a_z,A,'y',1.5);
fplot3(a_x'*T1,a_y'*T1,a_z'*T1,interv,'r','LineWidth',3);
title(['\textbf{Accel, long=',num2str(long_bs,4),'  $u_a$}'])

%% Other stuff

increm=0.001;
t_min=min(knots)+1.115;
t_max=max(knots)-1.4;

points_bs = bspline_deboor(p+1,knots,cPoints,t_min:increm:t_max);%0:increm:1


points_x=points_bs(1,:);
points_y=points_bs(2,:);
points_z=points_bs(3,:);

points_v_x=deltaT*diff(points_x)/(increm); 
points_v_y=deltaT*diff(points_y)/(increm);
points_v_z=deltaT*diff(points_z)/(increm);

points_a_x=deltaT*diff(points_v_x)/(increm);
points_a_y=deltaT*diff(points_v_y)/(increm);
points_a_z=deltaT*diff(points_v_z)/(increm);

length_pos=1;
length_vel=10;
length_accel=20;
view_pos=[260.4,-10.4];  %For MINVO paper, use [80,30];
view_vel=[19.73,50.49];  %For MINVO paper, use [136,30]
view_accel=[10.01,-14.19]; %For MINVO paper, use [136,30];

subplot(3,3,1);hold on;axis off; view(view_pos)
plot3(points_x,points_y,points_z,'-b'); plot_arrow_axes(length_pos); 
subplot(3,3,2);hold on;axis off; view(view_vel)
plot3(points_v_x,points_v_y,points_v_z,'-b'); plot_arrow_axes(length_vel)
subplot(3,3,3);hold on;axis off; view(view_accel)
plot3(points_a_x,points_a_y,points_a_z,'-b'); plot_arrow_axes(length_accel)


subplot(3,3,4);hold on;axis off; view(view_pos)
plot3(points_x,points_y,points_z,'-b'); plot_arrow_axes(length_pos)
subplot(3,3,5);hold on;axis off; view(view_vel)
plot3(points_v_x,points_v_y,points_v_z,'-b'); plot_arrow_axes(length_vel)
subplot(3,3,6);hold on;axis off; view(view_accel)
plot3(points_a_x,points_a_y,points_a_z,'-b'); plot_arrow_axes(length_accel)


subplot(3,3,7);hold on;axis off; view(view_pos)
plot3(points_x,points_y,points_z,'-b'); plot_arrow_axes(length_pos)
subplot(3,3,8);hold on;axis off; view(view_vel)
plot3(points_v_x,points_v_y,points_v_z,'-b'); plot_arrow_axes(length_vel)
subplot(3,3,9);hold on;axis off; view(view_accel)
plot3(points_a_x,points_a_y,points_a_z,'-b'); plot_arrow_axes(length_accel)


%%Export

%exportAsSvg(gcf,'comparisonBsBezierMinvoPosVelAccel_matlab_mader_paper')
%%

function plot_arrow_axes(length)

arrow3d([0 0 0],[0 0 length],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 length 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[length 0 0],20,'cylinder',[0.2,0.1]);

end