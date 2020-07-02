close all; clear; clc;

addpath(genpath('./utils'));

set(0,'DefaultFigureWindowStyle','docked') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


syms t real
interv=[-1,1];

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';


%%RESULT for 2D for a given polynomial
figure; hold on;
set(gcf, 'Position',  [500, 500, 3000, 1000])
% subplot(1,2,1);hold on

Ab=getSolutionA(2,"m11");


% P=[ 0.9349   -0.6   -0.9    0.5
%    -0.9349   -0.6    0.9    0.5];

novale=getSolutionA(3,"m11");
P=[novale(2,:); novale(3,:)];

Pb=P(:,2:end);
Pa=P(:,1);


V=Pb*inv(Ab);

% tmp=V';
% [k,av] = convhull(tmp);
% fill(tmp(k,1),tmp(k,2),'g','LineWidth',1)
% alpha 0.2
% pol_x=Pb(1,:)';
% pol_y=Pb(2,:)';
% fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);hold on;


% pol_x=[P(1,1) 0 0 0]';
% pol_y=[P(2,1) 0 0 0]';
% fplot(pol_x'*T3,pol_y'*T3,interv,'b','LineWidth',3);hold on;


Aa=inv([V;ones(1,3)])*[Pa; 0]

A=[Aa Ab];

% P-V*[Aa Ab]
% figure;
% fplot(A*T3, interv)


V=P*inv(getSolutionA(3,"m11"))
tmp=V';
[k,av] = convhull(tmp);
fill(tmp(k,1),tmp(k,2),'b','LineWidth',1)
alpha 0.2
axis equal
pol_x=P(1,:)';
pol_y=P(2,:)';
fplot(pol_x'*T3,pol_y'*T3,interv,'g','LineWidth',3);hold on;





tmp=inv(getSolutionA(3,"m11"));
P*tmp(:,1:3)


P=[ 0.1   -0.6   -0.9    0.5
   -0.9349   -0.6    0.9    0.5];

V=P*tmp(:,2:4)
tmp=V';
[k,av] = convhull(tmp);
fill(tmp(k,1),tmp(k,2),'r','LineWidth',1)
alpha 0.2
axis equal
pol_x=P(1,:)';
pol_y=P(2,:)';
fplot(pol_x'*T3,pol_y'*T3,interv,'r','LineWidth',3);hold on;

%%
figure; hold on;
clc;
A=getSolutionA(3,"m11");
fplot3(A(1,:)*T3,A(2,:)*T3,A(3,:)*T3,interv,'r','LineWidth',3);hold on;
d_x=polyder(A(1,:))
d_y=polyder(A(2,:))
d_z=polyder(A(3,:))
fplot3(d_x*T2,d_y*T2,d_z*T2,interv,'r','LineWidth',3);hold on;