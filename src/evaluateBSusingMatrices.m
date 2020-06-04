%Jesus Tordesillas Torres, June 2020

clc; clear; close all;
addpath(genpath('./utils'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


p = 3;
n = 4; %4 for degree 3
m=n+p+1;

global knots
knots = [0   0   0   0   0.1  0.2  0.3  0.4    0.5   0.5   0.5   0.5];  % knot vector

cPoints=[ 1         0.5     1.0    3.0   5.0    6.0     3           4;
          -1          0.3   1.0    1.0   2.0    0.0      1           1;
          -1          -1    2.0    1.0   0.0    0.0      0           1];
            
increm=0.001;

figure; hold on;

pos = bspline_deboor(n,knots,cPoints,knots(4):increm:knots(5) );%0:increm:1
plot3(pos(1,:),pos(2,:),pos(3,:),'r')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(5):increm:knots(6) );%0:increm:1
plot3(pos(1,:),pos(2,:),pos(3,:),'b')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(6):increm:knots(7) );%0:increm:1
plot3(pos(1,:),pos(2,:),pos(3,:),'g')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(7):increm:knots(8) );%0:increm:1
plot3(pos(1,:),pos(2,:),pos(3,:),'k')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(8):increm:knots(9) );%0:increm:1
plot3(pos(1,:),pos(2,:),pos(3,:),'y')%,'-b'

Abs=computeMatrixForBSpline(3,"01");

%%
syms u;

%First interval
last4cps= cPoints(:,1:4);
curve=last4cps*computeAforIntervalThatStartsAt_knots(4)*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'r--','LineWidth',3);

%Second interval
last4cps= cPoints(:,2:5);
curve=last4cps*computeAforIntervalThatStartsAt_knots(5)*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'b--','LineWidth',3);

%Third interval
last4cps= cPoints(:,3:6);
curve=last4cps*Abs*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'g--','LineWidth',3);

%Fourth interval
last4cps= cPoints(:,4:7);
curve=last4cps*computeAforIntervalThatStartsAt_knots(7)*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'k--','LineWidth',3);

%Fifth interval
last4cps= cPoints(:,5:8);
curve=last4cps*computeAforIntervalThatStartsAt_knots(8)*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'y--','LineWidth',3);

%%
clc
A_first_12=12*computeAforIntervalThatStartsAt_knots(4)
A_second_12=12*computeAforIntervalThatStartsAt_knots(5)

6*computeAforIntervalThatStartsAt_knots(6) %Same as 6*computeMatrixForBSpline(3,"01")

A_secondlast_12=12*computeAforIntervalThatStartsAt_knots(7) 
A_last_12=12*computeAforIntervalThatStartsAt_knots(8) %Same as 6*computeMatrixForBSpline(3,"01")


function A=computeAforIntervalThatStartsAt_knots(index_t_start)

%Following the notation from
%https://link.springer.com/article/10.1007/s003710050206 
%("General matrix representations for B-splines"
% See bottom box page 181

global knots
% i=3;

j=index_t_start;

ti=knots(j); tiP1=knots(j+1); tiP2=knots(j+2); tiP3=knots(j+3); tiM1=knots(j-1); tiM2=knots(j-2);

m00=((tiP1-ti)^2)/((tiP1-tiM1)*(tiP1-tiM2));
m02=((ti-tiM1)^2)/((tiP2 -tiM1)*(tiP1 -tiM1));
m01=1-m00-m02;
m03=0;

m10=-3*m00; 
m12 =3*(tiP1-ti)*(ti-tiM1)/((tiP2-tiM1)*(tiP1 -tiM1));
m11 =3*m00-m12;
m13=0;

m20 = 3*m00; 
m22 = 3*((tiP1-ti)^2)/((tiP2-tiM1)*(tiP1 -tiM1));
m21 = -3*m00-m22;
m23= 0;


m30 = -m00;
m33= ((tiP1-ti)^2)/((tiP3-ti)*(tiP2-ti));
m32 = -m22/3 - m33 - ((tiP1-ti)^2)/((tiP2-ti)*(tiP2-tiM1));
m31 = m00-m32-m33;
M=[m00 m01 m02 m03;
   m10 m11 m12 m13;
   m20 m21 m22 m23;
   m30 m31 m32 m33];

%And now change it to the convention I use

A=[M(4,:)' M(3,:)' M(2,:)' M(1,:)'];


end
