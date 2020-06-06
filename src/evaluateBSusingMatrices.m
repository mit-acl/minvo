%Jesus Tordesillas Torres, June 2020

clc; clear; close all;
addpath(genpath('./utils'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%%Degree=2
n = 3; %3 for degree 2


n_int_knots=4
deltaT=1/(n_int_knots+1);
interm=deltaT*(1:n_int_knots);
knots = [0   0   0   interm   max(interm)+deltaT max(interm)+deltaT max(interm)+deltaT];  % knot vector

cPoints=[ 1         0.5     1.0    3.0        3           2   5;
          -1          0.3   1.0    1.0           4           1  2];
            
increm=0.001;

figure; hold on;

pos = bspline_deboor(n,knots,cPoints,knots(3):increm:knots(4) );%0:increm:1
plot(pos(1,:),pos(2,:),'k')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(4):increm:knots(5) );%0:increm:1
plot(pos(1,:),pos(2,:),'r')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(5):increm:knots(6) );%0:increm:1
plot(pos(1,:),pos(2,:),'b')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(6):increm:knots(7) );%0:increm:1
plot(pos(1,:),pos(2,:),'g')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(7):increm:knots(8) );%0:increm:1
plot(pos(1,:),pos(2,:),'k')%,'-b'
pos = bspline_deboor(n,knots,cPoints,knots(8):increm:knots(9) );%0:increm:1
plot(pos(1,:),pos(2,:),'y')%,'-b'


syms u;


deg=2;
%First interval
last3cps= cPoints(:,1:3);
curve=last3cps*computeMatrixForClampedUniformBSpline(deg,0,"01")*[u*u u 1]';
fplot(curve(1),curve(2), [0,1],'k--','LineWidth',3);

disp("====")
%Second interval
last3cps= cPoints(:,2:4);
curve=last3cps*computeMatrixForClampedUniformBSpline(deg,1,"01")*[u*u u 1]';
fplot(curve(1),curve(2), [0,1],'r--','LineWidth',3);


disp("====")
%Third interval
last3cps= cPoints(:,3:5);
curve=last3cps*computeMatrixForClampedUniformBSpline(deg,2,"01")*[u*u u 1]';
fplot(curve(1),curve(2), [0,1],'b--','LineWidth',3);

%Fourth interval
last3cps= cPoints(:,4:6);
curve=last3cps*computeMatrixForClampedUniformBSpline(deg,-2,"01")*[u*u u 1]';
fplot(curve(1),curve(2), [0,1],'g--','LineWidth',3);

%Fifth interval
last3cps= cPoints(:,5:7);
curve=last3cps*computeMatrixForClampedUniformBSpline(deg,-1,"01")*[u*u u 1]';
fplot(curve(1),curve(2), [0,1],'y--','LineWidth',3);
%%
clc
2*computeMatrixForAnyBSpline(deg,3, knots,"01")
2*computeMatrixForAnyBSpline(deg,4, knots,"01")
2*computeMatrixForAnyBSpline(deg,5, knots,"01")
2*computeMatrixForAnyBSpline(deg,6, knots,"01")
2*computeMatrixForAnyBSpline(deg,7, knots,"01")

%%
clc
knots=[0 1 2 3 4 5]
deg=2
computeMatrixForAnyBSpline(deg,3, knots,"01");
%% Degree=3


n = 4; %4 for degree 3

knots = [0   0   0   0   0.1  0.2  0.3  0.4    0.5   0.5   0.5   0.5];  % knot vector

cPoints=[ 1         0.5     1.0    3.0   5.0    6.0     3           4;
          -1          0.3   1.0    1.0   2.0    0.0      1           1;
          -1          -1    2.0    1.0   0.0    0.0      0           1];
            
increm=0.001;

figure; hold on; view(90,90);
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

syms u;

deg=3;
%First interval
last4cps= cPoints(:,1:4);
curve=last4cps*computeMatrixForClampedUniformBSpline(deg,0,"01")*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'r--','LineWidth',3);

%Second interval
last4cps= cPoints(:,2:5);
curve=last4cps*computeMatrixForClampedUniformBSpline(deg,1,"01")*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'b--','LineWidth',3);

%Third interval
last4cps= cPoints(:,3:6);
curve=last4cps*Abs*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'g--','LineWidth',3);

%Fourth interval
last4cps= cPoints(:,4:7);
curve=last4cps*computeMatrixForClampedUniformBSpline(deg,-2,"01")*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'k--','LineWidth',3);

%Fifth interval
last4cps= cPoints(:,5:8);
curve=last4cps*computeMatrixForClampedUniformBSpline(deg,-1,"01")*[u*u*u u*u u 1]';
fplot3(curve(1),curve(2),curve(3), [0,1],'y--','LineWidth',3);

view(55,19);
%%
clc
A_first_12=12*computeMatrixForAnyBSpline(deg, 4, knots,"01")
A_second_12=12*computeMatrixForAnyBSpline(deg, 5, knots,"01")

6*computeMatrixForAnyBSpline(deg, 6, knots,"01") %Same as 6*computeMatrixForBSpline(3,"01")

A_secondlast_12=12*computeMatrixForAnyBSpline(deg, 7, knots,"01") 
A_last_12=12*computeMatrixForAnyBSpline(deg, 8, knots,"01") %Same as 6*computeMatrixForBSpline(3,"01")



