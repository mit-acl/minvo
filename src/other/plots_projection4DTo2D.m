%% Non-rational curves (projection)
close all; clear; clc;

addpath(genpath('./../utils'));
addpath(genpath('./../solutions'));

interv=[-1,1];

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');set(0,'defaultfigurecolor',[1 1 1])

syms t real

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';
V=[zeros(4,1) eye(4)]; %Standard simplex

linewidth=2;

A=getA_MV(4,interv);
% P=V*A;

T4=getT(4,t);

px=A(2,:)*T4;
py=A(3,:)*T4;
pz=A(4,:)*T4;
pw=A(5,:)*T4;
assume(t,'real')
figure;
subplot(2,2,1);hold on;

% curve_projected=(1/(px+pz))*[px pz]';

curve_projected=(1/(1-pw-py))*[px pz]';

patch([0 0 1],[0 1 0],'green','FaceAlpha',.3);
fplot(curve_projected(1),curve_projected(2),interv,'LineWidth',linewidth); %Note that we are using the standard simplex


A=getA_MV(3,interv);

px=A(2,:)*T3;
py=A(3,:)*T3;
pz=A(4,:)*T3;

curve_projected=(-1/(py-1))*[px pz]';

fplot(curve_projected(1),curve_projected(2),interv,'LineWidth',linewidth); %Note that we are using the standard simplex
