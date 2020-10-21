% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

close all; clear; clc;

addpath(genpath('./utils'));
addpath(genpath('./solutions'));

interv=[-1,1];

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%Let us change now the usual grey background of the matlab figures to white
%See https://www.mathworks.com/matlabcentral/answers/96816-how-do-i-change-the-default-background-color-of-all-figure-objects-created-in-matlab
set(0,'defaultfigurecolor',[1 1 1])

syms t real

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';

%% Print the ratios of the determinants:
disp('abs( det(A_MV)/det(A_Be) )')

for i=1:7
    value=abs(det(getA_MV(i,interv))/det(getA_Be(i,interv)));
    vpa(value,4)
end

disp('abs( det(A_MV)/det(A_BS) )')

for i=1:7
    value=abs(det(getA_MV(i,interv))/det(getA_BS(i,interv)));
    vpa(value,4)
end

%% Print all the roots

for i=1:7
    fprintf("n=%f\n",i)
    [A rootsA]=getA_MV(i,interv);
    matrix_with_roots=[];
    for j=1:size(rootsA,2)
        tmp=rootsA{j};
        if(length(tmp)<(i+1)/2)
            tmp=[tmp nan];
        end
        matrix_with_roots=[matrix_with_roots; tmp];
    end
    matrix_with_roots 
    latex(vpa(sym(matrix_with_roots),4));
end

%% Plot all the polynomials of the basis

name_figure='imgs/plots_basis';



figure;
n_rows=8;
n_cols=4;

subplot(n_rows,n_cols,1); hold on
set(gcf, 'Position',  [500, 500, 2000, 2500])



font_size_title=12;

for degree=1:4
   
   T=[];
    for i=0:(degree)
       T=[t^i ;T];
    end


   subplot(n_rows,n_cols,degree);
   fplot(getA_MV(degree,interv)*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{MINVO, n=',num2str(degree),'}'),'FontSize',font_size_title )
   box on
   
   subplot(n_rows,n_cols,n_cols+degree);
   fplot(getA_Be(degree,interv)*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{Bernstein, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
      
   subplot(n_rows,n_cols,2*n_cols+degree);
   fplot(getA_BS(degree,interv)*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{B-Spline, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,3*n_cols+degree);
   fplot(lagrangePoly(linspace(min(interv),max(interv),degree+1))*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{Lagrange, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
  
   
end


for degree=5:7 %TODO: Add 8
   
   T=[];
    for i=0:(degree)
       T=[t^i ;T];
    end

   subplot(n_rows,n_cols,4*n_cols+1+(degree-5));
   fplot(getA_MV(degree,interv)*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{MINVO, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,5*n_cols+1+(degree-5));
   fplot(getA_Be(degree,interv)*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{Bernstein, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   
   subplot(n_rows,n_cols,6*n_cols+1+(degree-5));
   fplot(getA_BS(degree,interv)*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{B-Spline, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,7*n_cols+1+(degree-5));
   fplot(lagrangePoly(linspace(min(interv),max(interv),degree+1))*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{Lagrange, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
  
      
end


% set(gca, 'Position',[0.7813, 0.1100, 0.0371, 0.8150]);

sp_hand1=subplot(n_rows,n_cols,[20, 24, 28])
plot(0,0,  0,0,  0,0,  0,0,  0,0,  0,0, 0,0 ,0,0,  0,0)
axis off
lgd=legend('$\lambda_0(t)$','$\lambda_1(t)$','$\lambda_2(t)$','$\lambda_3(t)$','$\lambda_4(t)$','$\lambda_5(t)$','$\lambda_6(t)$','$\lambda_7(t)$')
lgd.FontSize = 17;

lgd.Position = [0.9,0.9,1,0.3].*lgd.Position;

% pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
% new_pos1 = [1,1,0.1,1].*pos1 %smaller width
% set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


%exportAsPdf(gcf,name_figure)

%% RESULT for 2D for a given polynomial
figure; hold on;
set(gcf, 'Position',  [500, 500, 3000, 1000])
% subplot(1,2,1);hold on

A=getA_MV(2,interv);

v1=[0.1, 0.9];
v2=[0.4  1.0];
v3=[0,  0.35]; 
%TODO: use P=VA below, instead of component by component

vx=[v1(1)  v2(1)  v3(1)]';
vy=[v1(2)  v2(2)  v3(2)]';

V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'g','LineWidth',1)
pol_x=A'*vx;
pol_y=A'*vy;


fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',2);hold on;


%Bernstein
pol_x=pol_x+[0 0 0.5]';
A=getA_Be(2,interv);
vx=inv(A')*pol_x;
vy=inv(A')*pol_y;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'b','LineWidth',1)
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',2);


%B-Spline
pol_x=pol_x+[0 0 1.5]';
A=getA_BS(2,interv);
vx=inv(A')*pol_x;
vy=inv(A')*pol_y;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'y','LineWidth',1)
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',2);

alpha 0.2
axis equal
xlim([-2 5.5])
ylim([-1 2.0])

% exportAsSvg(gcf,'comparison2d_poly_given_matlab')

%% RESULT for 3D for a given polynomial
figure;
subplot(1,3,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

view1=30;
view2=30;


vx=[ 0    0.7  0.1 0.5]';
vy=[1.1   0.4  0.1 1.3]';
vz=[0.8   1  0 0]';

V=[vx'; vy'; vz'];


A=getA_MV(3,interv);
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

P=[pol_x'; pol_y'; pol_z']

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.017);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2); axis equal;

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

%%%% Bernstein
subplot(1,3,2); hold on; 
A=getA_Be(3,interv);

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2); axis equal;

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);


%%%% BSpline
subplot(1,3,3); hold on; 
A=getA_BS(3,interv);

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.2);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2); axis equal;

% a1=subplot(1,3,1);
% a2=subplot(1,3,2);
% a2=subplot(1,3,2);
% allYLim = get([a1 a2], {'YLim'});
% allYLim = cat(2, allYLim{:});
% set([a1 a2], 'YLim', [min(allYLim), max(allYLim)]);

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

%exportAsSvg(gcf,'imgs/comparison3d_poly_given_matlab')

%% RESULT for 2D for a given simplex
figure; hold on;
set(gcf, 'Position',  [500, 500, 3000, 1000])
subplot(1,3,1);hold on

A=getA_MV(2,interv);


v1=[0.5,  0.0];
v3=[-0.5,  0.0];
v2=[0.0,  3/4]; 


vx=[v1(1)  v2(1)  v3(1)]';
vy=[v1(2)  v2(2)  v3(2)]';

V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'g','LineWidth',1)
pol_x=A'*vx;
pol_y=A'*vy;

% plot(0,0,'*')
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);
alpha 0.2;xlim([-0.7,0.7]);ylim([-0.1,1.7]);axis equal;

%%%% Bezier
subplot(1,3,2);hold on
A=getA_Be(2,interv);
pol_x=A'*vx;
pol_y=A'*vy;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'b','LineWidth',1)
ylim([0 2.0])
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);
 alpha 0.2;xlim([-0.7,0.7]);ylim([-0.1,1.7]);axis equal;
%%%% BSpline
subplot(1,3,3);hold on
A=getA_BS(2,interv);
pol_x=A'*vx;
pol_y=A'*vy;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'y','LineWidth',1)
ylim([0 2.0])
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);
alpha 0.2;xlim([-0.7,0.7]);ylim([-0.1,1.7]);axis equal;

% exportAsSvg(gcf,'comparison2d_simplex_given_matlab')


%% GEOMETRIC INTERPRETATION OF LAMBDA_I
figure;
subplot(1,2,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

A=getA_MV(3,interv);
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

P=[pol_x'; pol_y'; pol_z']

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
 axis equal;
ylim([-1.5,1.5]);

arrow3d([0 0 0],[0 0 1.5],10,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1.5 0],10,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1.5 0 0],10,'cylinder',[0.2,0.1]);

% roots_poly=real(roots(A(1,:))); %real to avoid numerical approximations
% plotSphere(position, radius, color)

subplot(1,2,2);
hold on; 

v0=[0.1 0.5 0]';
v1=[1 0.2 0.5]';
v2=[0.8 1 0.4]';
v3=[0.1 0.7 1.0]';

V=[v0 v1 v2 v3];
vx=V(1,:)';
vy=V(2,:)';
vz=V(3,:)';

pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(100, 40)
axis equal;
ylim([-1.5,1.5]); 

% a1=subplot(1,2,1);
% a2=subplot(1,2,2);
% allYLim = get([a1 a2], {'YLim'});
% allYLim = cat(2, allYLim{:});
% set([a1 a2], 'YLim', [min(allYLim), max(allYLim)]);

arrow3d([0 0 0],[0 0 1],10,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],10,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],10,'cylinder',[0.2,0.1]);

%exportAsSvg(gcf,'imgs/geom_meaning_lambdai')

%% Result for 3D for a given simplex
A=getA_MV(3,interv);

figure;
subplot(2,2,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 2000])

v0=[1.1    -1/sqrt(3)   0]';
v1=[-0.5   -1/sqrt(3)   0.3]';
v2=[0     2/sqrt(3)   0.8]';
v3=[0.3     0           4/sqrt(6)]';

V=[v0 v1 v2 v3];
vx=V(1,:)';vy=V(2,:)';vz=V(3,:)';
pol_x=A'*vx;pol_y=A'*vy;pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.05);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);

%%Plot volume Bezier
subplot(2,2,2);hold on
A=getA_Be(3,interv);
volumen_bezier=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.05);

fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);

%%Plot volume B-Spline
subplot(2,2,[3,4]);hold on
A=getA_BS(3,interv);
V=P*inv(A);   %P=VA

volumen_bs=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.3);
subplot(2,2,[3,4]);hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);


view1=80;view2=30;

subplot(2,2,1);
title(['\textbf{MINVO, Volume=',num2str(volumen_minvo,4),'  $u^3$}'])
view(view1, view2)
axis equal;
xlim([-1,2]);
ylim([-2,2]);
plotAxesArrows(1);axis off

subplot(2,2,2);
title(['\textbf{B\''ezier, Volume=',num2str(volumen_bezier,4),'  $u^3$}'])
view(view1, view2)
axis equal;plotAxesArrows(1);axis off

subplot(2,2,[3,4]);
title(['\textbf{B-Spline, Volume=',num2str(volumen_bs,4),'  $u^3$}'])
view(view1, view2)
axis equal;
arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
axis off
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

% exportAsSvg(gcf,'imgs/comparisonBSBeMV_matlab')

%% RESULT for 3D for a given simplex
figure;
subplot(1,3,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

view1=54;
view2=-4.5;

v0=[0.1 0.5 0]';
v1=[1 0.2 0.5]';
v2=[0.8 1 0.4]';
v3=[0.1 0.7 1.0]';

V=[v0 v1 v2 v3];
vx=V(1,:)';
vy=V(2,:)';
vz=V(3,:)';

A=getA_MV(3,interv);
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

P=[pol_x'; pol_y'; pol_z'];

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
 axis equal;
ylim([-1.5,1.5]);

poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])


%%%%% Bernstein
subplot(1,3,2); hold on; 
A=getA_Be(3,interv);
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
axis equal;
ylim([-1.5,1.5]); 

poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])



%%%%% BSpline
subplot(1,3,3); hold on; 
A=getA_BS(3,interv);
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
axis equal;
ylim([-1.5,1.5]); 

poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])


% WORKS:
%print(gcf,'imgs/comparison3d_simplex_given_matlab','-dpng','-r1000')
 
% DON'T WORK:
% exportAsPdf(gcf,'imgs/comparison_convex_hull')
% saveas(gcf,'imgs/comparison_convex_hull.eps')
% addpath('./utils/plot2svg/plot2svg')
% plot2svg("temperature_standard.svg");
% printeps(get(gcf,'Number'),'imgs/comparison_convex_hull')
% saveas(gcf,'imgs/comparison_convex_hull.png')

%% Plot the contact points, and check the centroidal condition (for the case n=3)
figure;hold on; view1=54; view2=-4.5;

V=[zeros(3,1) eye(3)];
vx=V(1,:)';   vy=V(2,:)';  vz=V(3,:)';
A=getA_MV(3,interv);  
P=V*A;
pol_x=P(1,:)';  pol_y=P(2,:)';  pol_z=P(3,:)';

sizes_spheres=0.015;

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',sizes_spheres);
view(view1, view2);  axis equal; ylim([-1.5,1.5]);

%This epsilons are simply to avoid artificial effects with MATLAB plotting
%(which happen when two surfaces are very close to each other)
tol_for_visual=0.003
tol_for_visual2=0.013;
pol_x_tmp=(pol_x+[0 0 0 tol_for_visual]')/(1+tol_for_visual2)
pol_y_tmp=(pol_y+[0 0 0 tol_for_visual]')/(1+tol_for_visual2)
pol_z_tmp=(pol_z+[0 0 0 tol_for_visual]')/(1+tol_for_visual2)

poly=[pol_x_tmp'*T3,pol_y_tmp'*T3,pol_z_tmp'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,0.3)%,'FaceColor','cyan'
fplot3(pol_x_tmp'*T3,pol_y_tmp'*T3,pol_z_tmp'*T3,interv,'r','LineWidth',3);
camlight; lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])

mean1=mean(V(:,1:3),2);
mean2=mean(V(:,2:4),2);
mean3=mean(V(:,[1,2,4]),2);
mean4=mean(V(:,[1,3,4]),2);

plotSphere(mean1,sizes_spheres,'g')
plotSphere(mean2,sizes_spheres,'g')
plotSphere(mean3,sizes_spheres,'g')
plotSphere(mean4,sizes_spheres,'g')

v0=V(:,1); v1=V(:,2); v2=V(:,3); v3=V(:,4);

all_roots=getAllRoots_MV(3,interv);
for roots_i=all_roots
    pt=double(P*subs(T3,t,roots_i));
    plotSphere(pt,sizes_spheres,'r')
end

ta=all_roots(end-1);
tb=all_roots(end-2);

Pm1=double(P*subs(T3,t,-1));
Pmtb=double(P*subs(T3,t,-tb));
Pmta=double(P*subs(T3,t,-ta));
Ptb=double(P*subs(T3,t,tb));
Pta=double(P*subs(T3,t,ta));
P1=double(P*subs(T3,t,1));

plotsegment(P1,Pmta,'c',1)
plotsegment(P1,Ptb,'c',1)
plotsegment(Pm1,Pta,'c',1)
plotsegment(Pm1,Pmtb,'c',1)

%Plot the legend (by hand)
s1=[0 0.5 2];
s2=[0 0.5 1.9];
s3=[0 0.5 1.8];
plotSphere(s1,0.02,'r')
plotSphere(s2,0.02,'g')
plotSphere(s3,0.02,[.98 .45 .02])
text(s1(1),s1(2)+0.05,s1(3), "Contact points");
text(s2(1),s2(2)+0.05,s2(3), "Centroid of each face");
text(s3(1),s3(2)+0.05,s3(3), "Control points");

Vsmall=[mean1 mean2 mean3 mean4];
[k1,av1] = convhull(Vsmall(1,:)',Vsmall(2,:)',Vsmall(3,:)');
trisurf(k1,Vsmall(1,:)',Vsmall(2,:)',Vsmall(3,:)','FaceColor','yellow','EdgeColor','k','FaceAlpha' ,1.0)

% print('-dpng','-r500',"centroid_facets")

%% RESULT for 3D for a given polynomial (with and without splitting)

num_of_intervals=5;

figure; last_figure=get(gcf,'Number');
set(gcf, 'Position',  [500, 500, 3000, 1000])
figure(last_figure+1); 
set(gcf, 'Position',  [500, 500, 3000, 1000])

figure(last_figure); 
subplot(1,4,1);hold on
view1=30; view2=30;
vx=[ 0    0.7  0.1 0.5]'; 
vy=[1.1   0.4  0.1 1.3]';
vz=[0.8   1  0 0]';
V=[vx'; vy'; vz'];
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.017)
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
plotAxesArrows(0.5)
view(view1, view2); axis equal;

figure(last_figure+1); 
subplot(1,4,1); hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
[volumen_minvo_splitted,num_vertexes]=plot_splitted_convex_hulls(P,A,interv,num_of_intervals,'g',0.017);
title(["volMVsplitted/volMV= ",num2str(volumen_minvo_splitted/volumen_minvo)," numVertexes= ",num2str(num_vertexes)]);
plotAxesArrows(0.5);
view(view1, view2); axis equal;


figure(last_figure);
subplot(1,4,2);hold on
view1=30; view2=30;
vx=[ 0.8    1.6  0.1 0.5]';
vy=[-2.3   -0.6  1.4 0.3]';
vz=[0.2   -1  0 0.7]';
V=[vx'; vy'; vz'];
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.037)
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
plotAxesArrows(1)
view(view1, view2); axis equal;

figure(last_figure+1); 
subplot(1,4,2); hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
[volumen_minvo_splitted,num_vertexes]=plot_splitted_convex_hulls(P,A,interv,num_of_intervals,'g',0.030)
title(["volMVsplitted/volMV= ",num2str(volumen_minvo_splitted/volumen_minvo)," numVertexes= ",num2str(num_vertexes)]);
plotAxesArrows(1);
view(view1, view2); axis equal;

figure(last_figure);
subplot(1,4,3);hold on
view1=73.4; view2=35.36;
vx=[ -1.2    0.2  2.3 0.1]';
vy=[0.3   -0.5  0.1 -0.5]';
vz=[0.7   1  -0.4 0.1]';
V=[vx'; vy'; vz'];
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.037)
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
plotAxesArrows(1)
view(view1, view2); axis equal;

figure(last_figure+1); 
subplot(1,4,3); hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
[volumen_minvo_splitted,num_vertexes]=plot_splitted_convex_hulls(P,A,interv,num_of_intervals,'g',0.025)
title(["volMVsplitted/volMV= ",num2str(volumen_minvo_splitted/volumen_minvo)," numVertexes= ",num2str(num_vertexes)]);
plotAxesArrows(1);
view(view1, view2); axis equal;

figure(last_figure);
subplot(1,4,4);hold on
view1=30; view2=30;
vx=[ -1.2    0.6  1.3 0.5]';
vy=[0.3   -0.7  -0.4  0.5]';
vz=[-1.2   1  -0.4 0.1]';
V=[vx'; vy'; vz'];
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_minvo=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.037)
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
plotAxesArrows(1)
view(view1, view2); axis equal;

figure(last_figure+1); 
subplot(1,4,4); hold on
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
[volumen_minvo_splitted,num_vertexes]=plot_splitted_convex_hulls(P,A,interv,num_of_intervals,'g',0.025)
title(["volMVsplitted/volMV= ",num2str(volumen_minvo_splitted/volumen_minvo)," numVertexes= ",num2str(num_vertexes)]);
plotAxesArrows(1);
view(view1, view2); axis equal;

% figure(last_figure); 
% exportAsSvg(gcf,'imgs/many_comparisons3d_poly_given_matlab')

figure(last_figure+1); 
% exportAsSvg(gcf,'imgs/splitted_many_comparisons3d_poly_given_matlab')

%% RESULT for 3D for a given simplex

figure;
subplot(1,4,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])
view1=135; view2=-10.85;
V=[   -0.6330   -0.2630    0.2512    0.5605;
   -0.8377    0.8588    0.5514   -0.0264;
   -0.1283   -0.1064   -0.3873    0.0170];
vx=V(1,:)'; vy=V(2,:)'; vz=V(3,:)';
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.02);
view(view1, view2); axis equal;
poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])
% 
subplot(1,4,2);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])
view1=289.97; view2=-4.4981;
V=[    0.1741   -0.0582   -0.6105   -0.5447;
   -0.5845   -0.5390   -0.5482   -0.1286;
   -0.3975    0.6886   -0.6586   -0.3778];
vx=V(1,:)'; vy=V(2,:)'; vz=V(3,:)';
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.010);
view(view1, view2); axis equal;
poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])

subplot(1,4,3);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])
view1=286.97; view2=0.29;
V=[    0.9234    0.9049    0.1111    0.5949;
    0.4302    0.9797    0.2581    0.2622;
    0.1848    0.4389    0.4087    0.6028];
vx=V(1,:)'; vy=V(2,:)'; vz=V(3,:)';
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.010);
view(view1, view2); axis equal;
poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])

subplot(1,4,4);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])
view1=286.97; view2=0.29;
V=[    0.7112    0.2967    0.5079    0.8010;
    0.2217    0.3188    0.0855    0.0292;
    0.1174    0.4242    0.2625    0.9289];
vx=V(1,:)'; vy=V(2,:)'; vz=V(3,:)';
A=getA_MV(3,interv);
pol_x=A'*vx; pol_y=A'*vy; pol_z=A'*vz;
P=[pol_x'; pol_y'; pol_z'];
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.006);
view(view1, view2); axis equal;
poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=min(interv):0.01:max(interv);
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])


%print(gcf,'imgs/many_comparisons3d_simplex_given_matlab','-dpng','-r1000')
%% Video

h=figure(7); %Select the figure you want
for i=1:size(h.Children)
    subplot(h.Children(i))
    axis off
     tmp=gca;
     if (size(findobj(tmp.Children,'Type','Light'))<1) %If still no light in the subplot
         camlight %create light
     end
    title("");
    axis vis3d
    OptionZ.FrameRate=30;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    %Uncomment next line to record and save the video%%%%%%%%%%%%%%%%%%%
    %CaptureFigVid([-20,10;-110,10;-190,10;-290,10;-380,10], ['./videos/comparison_given_curve_only_traj_',num2str(i)],OptionZ)
end

%[Does NOT work] To move all the subplots at the same time. But if used with vis3d,
%subplots overlaf
%Link = linkprop(h.Children, {'CameraUpVector', 'CameraPosition', 'CameraTarget'}) %, 'CameraTarget'
%setappdata(gcf, 'StoreTheLink', Link);


%% Plot the root distribution of the MINVO basis functions
figure; hold on;   set(gcf, 'Position',  [1726, 663, 635, 355])

color_MV=[123,218,104]/255; %[178,238,166]/255;

for n=1:7
    % subplot(7,1,n);  
    MV_points=getAllRoots_MV(n,interv);

    color=[0.4940, 0.1840, 0.5560];
    delta=0.05;
    y_pos=delta*(7-n+1);
    plot(MV_points,y_pos*ones(size(MV_points)),'-o','MarkerEdgeColor','k', 'MarkerFaceColor', color_MV,'Color',color_MV,'MarkerSize',9,'LineWidth',1.2);
    h = gca; h.YAxis.Visible = 'off';
    ylim([0,delta*(n+0.6)])
    text(-1.2,y_pos,(['\textbf{n=',num2str(n),'}']),'FontSize',11)

%     Uncomment the part below if you want to plot the LGL, LGR and LG points for comparison 
%     fplot(A*getT(n,t),interv,'k');
%     offset=0;
%     N=numel(MV_points)%2*n%+offset
%    %See page 46 of https://www.researchgate.net/publication/38002446_Advancement_and_analysis_of_Gauss_pseudospectral_transcription_for_optimal_control_problems
%     lgl=(1-t*t)*diff(legendreP(N-1,t),t);  %its roots are the LGL points
%     lgr=legendreP(N,t)-legendreP(N-1,t);  %its roots are the LGR points
%     lg=legendreP(N,t);%its roots are the LG points
% 
%     lgl_points=double(solve(lgl==0,t))
%     lgr_points=double(solve(lgr==0,t))
%     lg_points=double(solve(lg==0,t))
    
%     plot(lgl_points,0.4*ones(size(lgl_points)),'-o', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
%     plot(lgr_points,0.3*ones(size(lgr_points)),'-o', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
%     plot(lg_points,0.2*ones(size(lg_points)),'-o', 'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
%      title(strcat('\textbf{n=',num2str(n),'}'),'FontSize',10 , 'Position', [-1.1, 0.2, 0])
%     if(n==1)
%         legend({'LGL','LGL','LG','Minvo'})
%     end
     
end

title('\textbf{Roots of the MINVO basis functions}'); xlabel('t');
% exportAsPdf(gcf,'roots_distribution');


%% Embeddings with k=3, different n and m
figure; hold on; set(gcf, 'Position',  [500, 500, 3000, 2000])
A=getA_MV(3,interv);
%%Curve with n=3, with m=2 and k=3
P=[    0.5670   -0.7308   -0.1584    0.8993;
   -2.0953    0.3654    0.7380   -0.4995;
    0.7641    0.1827   -0.2898    0.3001];
P(2,:)=[0 0 0 1]-P(1,:)-2.0*P(3,:) %Plane x+y+2*z=1
subplot(1,3,1);hold on;
volumen_minvo=plot_convex_hull(P(1,:)', P(2,:)', P(3,:)',A,'g',0.03);
fplot3(P(1,:)*T3,P(2,:)*T3,P(3,:)*T3,interv,'r','LineWidth',3);
xlabel('x'); ylabel('y'); zlabel('z')
title(['\textbf{Pos, vol=',num2str(volumen_minvo,4),'  $u_p^3$}']); axis off;
plotAxesArrows(1.8);view(165,22)

%%Curve with n=2, with m=2 and k=3
P=0.35*[ -2.2849    1.7357    0.3745;
   -3.5932   -1.0549    1.3875;
    0.8552    0.4869    0.4401];
subplot(1,3,2); hold on;
A=getA_MV(2,interv);
area_minvo= plot_plane_convex_hull(P(1,:)', P(2,:)', P(3,:)', A, 'g', 0.04);
fplot3(P(1,:)*T2,P(2,:)*T2,P(3,:)*T2,interv,'r','LineWidth',3);
title(['\textbf{area=',num2str(area_minvo,4),'  $u_v^2$}']); axis off;
plotAxesArrows(1);view(-187,7)

%%Curve with n=1, with m=1 and k=3
subplot(1,3,3); hold on; 
P=0.05*[   -4.5697    1.7357;
   -7.1863   -1.0549;
    7.7104    7.4869]
A=getA_MV(1,interv);
long_minvo=plot_line_convex_hull(P(1,:)', P(2,:)', P(3,:)',A,'g',0.03);
fplot3(P(1,:)*T1,P(2,:)*T1,P(3,:)*T1,interv,'r','LineWidth',3);
title(['\textbf{long=',num2str(long_minvo,4),'  $u_a$}']); axis off;
plotAxesArrows(1);view(93,10)

% print('-dpng','-r500',"embedded_curves_matlab")

%% Non-rational curves (projection)

V=[zeros(3,1) eye(3)]; %Standard simplex

linewidth=2;

A=getA_MV(3,interv);
P=V*A;

% figure;
% vol=plot_convex_hull(P(1,:)',P(2,:)',P(3,:)',A,'g',0.0);
% fplot3(P(1,:)*T3,P(2,:)*T3,P(3,:)*T3,interv,'r','LineWidth',1);
% fplot3(curve_projected(1),curve_projected(2),curve_projected(3),interv,'--','LineWidth',2);

figure;
subplot(2,2,1);hold on;

v_proj=V(:,3);   V_rest=[V(:,1),V(:,2),V(:,4)];
curve_projected=projectPoint(P*T3,v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'w');
patch(V_rest(1,:),V_rest(3,:),'green','FaceAlpha',.3); 
fplot(curve_projected(1),curve_projected(3),interv,'r','LineWidth',linewidth); %Note that we are using the standard simplex
axis equal; axis off; title('$\pi_2$');

subplot(2,2,2);hold on;
v_proj=V(:,2);   V_rest=[V(:,1),V(:,3),V(:,4)];
curve_projected=projectPoint(P*T3,v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'w');
patch(V_rest(2,:),V_rest(3,:),'green','FaceAlpha',.3);
fplot(curve_projected(2),curve_projected(3),interv,'r','LineWidth',linewidth,'MeshDensity',300); %Note that we are using the standard simplex
axis equal; axis off;title('$\pi_1$');

subplot(2,2,3);hold on;
v_proj=V(:,4);   V_rest=[V(:,1),V(:,2),V(:,3)];
curve_projected=projectPoint(P*T3,v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'w');
patch(V_rest(1,:),V_rest(2,:),'green','FaceAlpha',.3); 
fplot(curve_projected(1),curve_projected(2),interv,'r','LineWidth',linewidth,'MeshDensity',300); %Note that we are using the standard simplex
axis equal; axis off; title('$\pi_3$');

subplot(2,2,4); hold on;
v_proj=V(:,1);   V_rest=[V(:,2),V(:,3),V(:,4)];
curve_projected=projectPoint(P*T3,v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'c');
V_rest1_camera_plane=projectPoint(V_rest(1,:)',v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'c');
V_rest2_camera_plane=projectPoint(V_rest(2,:)',v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'c');
V_rest3_camera_plane=projectPoint(V_rest(3,:)',v_proj,V_rest(:,1),V_rest(:,2),V_rest(:,3),'c');
tmp=[V_rest1_camera_plane  V_rest2_camera_plane V_rest3_camera_plane];
patch(tmp(1,:),tmp(2,:),'green','FaceAlpha',.3);
fplot(curve_projected(1),curve_projected(2),interv,'r','LineWidth',linewidth,'MeshDensity',300); 
axis equal; axis off;title('$\pi_0$');

% exportAsSvg(gcf,'projections_matlab');

%% CONE OF POSITIVE POLYNOMIALS FOR n=2

figure;  hold on;
n=0.05%Change to .005 for a finer resolution;
size_grid=1.2%4.7; %was 1.2
% create coordinates
[a,b,c] = meshgrid(-size_grid:n:size_grid,-size_grid:n:size_grid,-size_grid:n:size_grid);

increm=0.001;
interv=[-1,1];
tt=min(interv):increm:max(interv);
region3=1.0;

syms t real

for ti=tt
    ti
%     root(aa*t^2+bb*t+cc,t)
    region3=region3 & (a*ti.^2+b*ti+c)>=0;
end

p=patch(isosurface(a,b,c,region3,0.0));
set(p,'FaceColor',[255,169,107]/255,'EdgeColor','none','FaceAlpha',.9);
  xlabel('a');ylabel('b');zlabel('c');
camlight 
lighting gouraud
grid on

%%%
% pos_root=(-b+sqrt(b.*b-4*a.*c))./(2*a);
% neg_root=(-b-sqrt(b.*b-4*a.*c))./(2*a);
% region3=(( imag(pos_root)==0 &(pos_root<=max(interv)) & (pos_root>=min(interv)) ) | ...
%         (imag(neg_root)==0 &(neg_root<=max(interv)) & (neg_root>=min(interv))));
% region3=double(region3);
% region3(region3>0.1)=100;
% p=patch(isosurface(a,b,c,region3,0.0));
%%%

A2=getA_MV(2,interv);

arrow3d([0 0 0],A2(1,:),10,'cylinder',[0.4,0.4]);
arrow3d([0 0 0],A2(2,:),10,'cylinder',[0.2,0.45]);
arrow3d([0 0 0],A2(3,:),10,'cylinder',[0.4,0.45]);

V=[0;0;0];
V=[V A2(1,:)' A2(2,:)' A2(3,:)'  A2(1,:)'+A2(2,:)' A2(2,:)'+A2(3,:)' A2(1,:)'+A2(3,:)'  A2(1,:)'+A2(2,:)'+A2(3,:)'];

pts=[V(:,1),V(:,2)];plot3(pts(1,:),pts(2,:),pts(3,:));
pts=[V(:,1),V(:,3)];plot3(pts(1,:),pts(2,:),pts(3,:));
pts=[V(:,1),V(:,4)];plot3(pts(1,:),pts(2,:),pts(3,:));
pts=[V(:,5),V(:,8)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,6),V(:,8)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,7),V(:,8)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,2),V(:,7)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,2),V(:,5)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,3),V(:,5)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,3),V(:,6)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,4),V(:,6)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');
pts=[V(:,4),V(:,7)];plot3(pts(1,:),pts(2,:),pts(3,:),'k');

[k1,av1] = convhull(V(1,:),V(2,:),V(3,:));
trisurf(k1,V(1,:),V(2,:),V(3,:),'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.15)
set(gca,'DataAspectRatio',[4 6 1])
view([88,29])
axis equal; grid off; axis off;

arrow3D_v2([0,0,0] ,sum(A2), 'b', 0.75,0.02);
arrow3D_v2([0,0,0] ,[1,0,0], 'r', 0.75,0.02);
arrow3D_v2([0,0,0] ,[0,1,0], 'g', 0.75,0.02);

% print('-dpng','-r500',"cone_A2")

%This plots the boundary of the polynomials that are positive for ALL t
% hold on;
% f = @(aa,bb,cc)  bb^2 - 4*aa*cc;
% interval = [-1.5 1.5 -1.5 1.5 0 1.5];
% fimplicit3(f,interval)

%This plots the coeff of the polynomials that are (t-r1)^2
% syms r1 t real
% fplot3(sym(1),-2*r1,r1*r1,interv)



% A2_scaled=A2; 
% A2_scaled(:,1)=A2_scaled(:,1)./A2(:,3);
% A2_scaled(:,2)=A2_scaled(:,2)./A2(:,3);
% A2_scaled(:,3)=A2_scaled(:,3)./A2(:,3);
% arrow3d([0,0,0] ,scaled(1,:),10,'cylinder',[0.4,0.4]);
% arrow3d([0,0,0] ,scaled(2,:),10,'cylinder',[0.4,0.4]);
% arrow3d([0,0,0] ,scaled(3,:),10,'cylinder',[0.4,0.4]);
% 
% sp=4.0;
% [x y] = meshgrid(-sp:0.1:sp); % Generate x and y data
% z = 1.0+0.0*x; % Generate z data
% surf(x, y, z,'FaceAlpha',0.15) % Plot the surface

%% SURFACES!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start code from https://www.mathworks.com/matlabcentral/fileexchange/37876-construction-of-cubic-bezier-patch-and-surface
load('teapot'); %loading matrix S. The file teapot.mat is available at BezierPatchSurface/BezierPatchSurface
% % Matrix S stores all the control points of all the CUBIC patches of
% % teapot surface such that
% % S(:,:,:,k) control points of kth patch, where k=1..32
% % Size of S(:,:,:,k) is 4 x 4 x 3, i.e., 16 control points and each
% % control point has three values (x,y,z)

% % S(:,:,1,k): x-coordates of control points of kth patch as 4 x 4 matrix 
% % S(:,:,2,k): y-coordates of control points of kth patch as 4 x 4 matrix 
% % S(:,:,3,k): z-coordates of control points of kth patch as 4 x 4 matrix
% % ------------------------------------
[r c d np]=size(S);
% % np: number of patches
ni=40; %number of interpolated values between end control points
u=linspace(0,1,ni); v=u;  %uniform parameterization
% % Higher the value of ni smoother the surface but computationally
% % expensive
% % ------------------------------------
% % Cubic Bezier interpolation of control points of each patch

for k=1:np
    Q(:,:,:,k)=bezierpatchinterp(S(:,:,:,k),u,v); %interpolation of kth patch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End code from https://www.mathworks.com/matlabcentral/fileexchange/37876-construction-of-cubic-bezier-patch-and-surface


k_plot=[5,15,21];
teapot_fig=figure; hold on;
for k=1:np
    if(ismember(k,k_plot))
         surface(Q(:,:,1,k),Q(:,:,2,k),Q(:,:,3,k),'FaceColor','red','EdgeColor','none')
    else
         surface(Q(:,:,1,k),Q(:,:,2,k),Q(:,:,3,k),'FaceColor','green','EdgeColor','none')
    end
end

camlight; axis equal; axis off;
view(3); box;  view1=12; view2=18.09; view(view1,view2)

% print(teapot_fig,'-dpng','-r500',"teapot_matlab")

n=3; %They are cubic patches
m=3;

An_Be=getA_Be(n,interv);
Am_Be=getA_Be(m,interv);

An_MV=getA_MV(n,interv);
Am_MV=getA_MV(m,interv);


figure;set(gcf, 'Position',  [500, 500, 2000, 2500]);tiledlayout('flow')

for i_plot=1:length(k_plot)
    
k=k_plot(i_plot);    

SBe_x=S(:,:,1,k);%S(:,:,1,k) has the x-coordinates of all the Bezier control points of the k-th patch
SBe_y=S(:,:,2,k);%S(:,:,2,k) has the y-coordinates of all the Bezier control points of the k-th patch
SBe_z=S(:,:,3,k);%S(:,:,3,k) has the z-coordinates of all the Bezier control points of the k-th patch

SMV_x=inv(An_MV)'*An_Be'*SBe_x*Am_Be*inv(Am_MV); %x-coordinates of all the MINVO control points of the k-th patch
SMV_y=inv(An_MV)'*An_Be'*SBe_y*Am_Be*inv(Am_MV); %...
SMV_z=inv(An_MV)'*An_Be'*SBe_z*Am_Be*inv(Am_MV); %...

syms u v real
Tn=[];
for i=0:(n)
     Tn=[u^i ;Tn];
end

Tm=[];
for i=0:(m)
     Tm=[v^i ;Tm];
end

%MINVO
nexttile;hold on;axis equal; view(view1,view2)
%See https://cse.taylor.edu/~btoll/s99/424/res/ucdavis/CAGDNotes/Matrix-Cubic-Bezier-Patch/Matrix-Cubic-Bezier-Patch.html
surface=sym(zeros(3,1));
surface(1,:)=Tn'*An_MV'*SMV_x*Am_MV*Tm; 
surface(2,:)=Tn'*An_MV'*SMV_y*Am_MV*Tm; 
surface(3,:)=Tn'*An_MV'*SMV_z*Am_MV*Tm; 
fsurf(surface(1), surface(2), surface(3), [min(interv), max(interv), min(interv), max(interv)],'r'); title("")

[k1,vol_MV] = convhull(SMV_x(:),SMV_y(:),SMV_z(:));
trisurf(k1,SMV_x(:),SMV_y(:),SMV_z(:),'FaceColor','g','FaceAlpha',0.2 )
camlight;material shiny ; axis off;

%Bernstein (\equiv Bezier)
%MINVO
nexttile;hold on;axis equal; view(view1,view2)
%See https://cse.taylor.edu/~btoll/s99/424/res/ucdavis/CAGDNotes/Matrix-Cubic-Bezier-Patch/Matrix-Cubic-Bezier-Patch.html
surface=sym(zeros(3,1));
surface(1,:)=Tn'*An_Be'*SBe_x*Am_Be*Tm; 
surface(2,:)=Tn'*An_Be'*SBe_y*Am_Be*Tm; 
surface(3,:)=Tn'*An_Be'*SBe_z*Am_Be*Tm; 
fsurf(surface(1), surface(2), surface(3), [min(interv), max(interv), min(interv), max(interv)],'r'); title("");   %Should be the same curve as before

[k1,vol_Be] = convhull(SBe_x(:),SBe_y(:),SBe_z(:));
trisurf(k1,SBe_x(:),SBe_y(:),SBe_z(:),'FaceColor','b','FaceAlpha',0.2 )
%Alternative form: see https://en.wikipedia.org/wiki/B%C3%A9zier_surface#Equation
% surface=Tn*An*K
% surface=zeros(3,1);
% total=0;
% for i=1:n+1
%     for j=1:m+1
%        total=total+1;
%        surface= surface+ An(i,:)*Tn*Am(j,:)*Tm*KMV(:,total);
%     end
% end


camlight; material shiny; axis off;

sprintf('vol_Be/vol_MV= %f',vol_Be/vol_MV)


axis equal

end

%% Different degrees of curves in 2D and 3D

figure; hold on; tiledlayout('flow');

%BOXPLOTS FOR 2D CURVES, WITH DIFFERENT n
a = -1;
b = 1;

for deg=2:7
    ratios=[];
    A_MV=getA_MV(deg,interv); A_MV_inv=A_MV^(-1);
    A_Be=getA_Be(deg,interv); A_Be_inv=A_Be^(-1);
    
    tt=linspace(min(interv),max(interv),deg+1);
    
    for (i=1:10000)
        i
      
    P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
        V_MV=P*A_MV_inv;
        V_Be=P*A_Be_inv;
        
       [k,aMV] = convhull(V_MV');       
       [k,aBe] = convhull(V_Be'); 
       ratios=[ratios aBe/aMV];
    end
    nexttile; hold on;
    
%     histogram(ratios,300,'Normalization','probability'); 
    boxplot(ratios,'symbol',''); yl = ylim; ylim([0.0,max(yl)]);
    if(deg==2)
        ylim([0,2*max(yl)]);
    end
    if(deg==7)
        ylim([-3.0,max(yl)]);
    end   
    yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); ylabel('$r$','FontSize', 12);  %'LineSpec','-.'
    title(['\textbf{n=',num2str(deg),'}'],'FontSize', 17);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    set(gcf,'Position',[1,1,1901,225]);
end
%RANDOM CASES FOR 2D CURVES, WITH DIFFERENT n
for (deg=2:7)
    nexttile; hold on;
    P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
    A_MV=getA_MV(deg,interv); 
    A_Be=getA_Be(deg,interv); 
        V_MV=P*A_MV^(-1);
        V_Be=P*A_Be^(-1);
        
       [k,aMV] = convhull(V_MV');       
       [k,aBe] = convhull(V_Be'); 
       
    [k,aMV] = convhull(V_MV'); conv_minvo=plot(V_MV(1,k),V_MV(2,k),'Color',[73 204 73]/255,'LineWidth',2);
    [k,aBe] = convhull(V_Be'); conv_Be=plot(V_Be(1,k),V_Be(2,k),'Color',[98 98 200]/255,'LineWidth',2);
    fplot(P(1,:)*getT(deg,t),P(2,:)*getT(deg,t),interv,'r','LineWidth',2); 
    xlabel('$x(t)$'); ylabel('$y(t)$'); title(['\textbf{n=',num2str(deg),'}',', $r=$',num2str(aBe/aMV,3)],'FontSize', 17)
end
hL = legend([conv_minvo,conv_Be],{'MINVO',"B\'{ezier}"},'FontSize', 12);
set(gcf,'Position',[166,880,2329,336])
% exportAsPdf(gcf,'diff_degree2D');
% savefig('diff_degree2D');

%BOXPLOTS FOR 3D CURVES, WITH DIFFERENT n

figure; hold on; tiledlayout('flow');
for deg=3:7
    ratios=[];
    A_MV=getA_MV(deg,interv); A_MV_inv=A_MV^(-1);
    A_Be=getA_Be(deg,interv); A_Be_inv=A_Be^(-1);
    tt=linspace(min(interv),max(interv),deg+1);
    for (i=1:10000)
        i
      
    P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
        V_MV=P*A_MV_inv;
        V_Be=P*A_Be_inv;
        
       [k,vol_MV] = convhull(V_MV(1,:),V_MV(2,:),V_MV(3,:));       
       [k,vol_Be] = convhull(V_Be(1,:),V_Be(2,:),V_Be(3,:)); 
       ratios=[ratios vol_Be/vol_MV];
    end
    nexttile; hold on;
    
%     histogram(ratios,300,'Normalization','probability'); 
    boxplot(ratios,'symbol',''); yl = ylim; ylim([0.0,max(yl)]);
    if(deg==3)
        ylim([0,2*max(yl)]);
    end
    if(deg==7)
        ylim([-3.0,max(yl)]);
    end   
    yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); ylabel('$r$','FontSize', 12);  %'LineSpec','-.'
    title(['\textbf{n=',num2str(deg),'}'],'FontSize', 17);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
end
set(gcf,'Position',[166,880,2329,336])
% exportAsPdf(gcf,'diff_degree3D_boxplots');
% savefig('diff_degree3D_boxplots');

%RANDOM CASES FOR 3D CURVES, WITH DIFFERENT n
figure; hold on; 
j=1;
size_arrows=0.2;
hlinks=[];
for (deg=3:7)
    ax1=subplot(2,5,j); j=j+1;
     hold on;
     a=0.0; b=0.2;

    P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
    A_MV=getA_MV(deg,interv); 
    A_Be=getA_Be(deg,interv); 
        V_MV=P*A_MV^(-1);
        V_Be=P*A_Be^(-1);

    [k,vol_MV] = convhull(V_MV(1,:),V_MV(2,:),V_MV(3,:)); trisurf(k,V_MV(1,:),V_MV(2,:),V_MV(3,:),'FaceColor','g','FaceAlpha',0.2);

    fplot3(P(1,:)*getT(deg,t),P(2,:)*getT(deg,t),P(3,:)*getT(deg,t),interv,'r','LineWidth',2); 
    camlight; axis equal; axis off;plotAxesArrows(size_arrows); lighting phong;

    ax2=subplot(2,5,j); j=j+1;
    hold on;
    [k,vol_Be] =  convhull(V_Be(1,:),V_Be(2,:),V_Be(3,:)); trisurf(k,V_Be(1,:),V_Be(2,:),V_Be(3,:),'FaceColor','b','FaceAlpha',0.2);
    fplot3(P(1,:)*getT(deg,t),P(2,:)*getT(deg,t),P(3,:)*getT(deg,t),interv,'r','LineWidth',2); 
    xlabel('$x(t)$'); ylabel('$y(t)$'); title(['\textbf{n=',num2str(deg),'}',', $r=$',num2str(vol_Be/vol_MV,3)]);
    camlight; axis equal; axis off;plotAxesArrows(size_arrows); lighting phong;

    hlinks = [hlinks linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'})];
end

% print('-dpng','-r500',"diff_degree3D_matlab")

%% Comparison with the area/volume of the convex hull of the curve

%%%%%%%%%%%%%%%%%%%%Case 2D (k=m=2)
means_MV=[]; means_Be=[];
stds_MV=[]; stds_Be=[];

all_degs=2:7;
for deg=all_degs
    
    ratios_MV=[];
    ratios_Be=[];
    
    A_MV=getA_MV(deg,interv); A_MV_inv=A_MV^(-1);
    A_Be=getA_Be(deg,interv); A_Be_inv=A_Be^(-1);
    tt=linspace(min(interv),max(interv),deg+1);
    for (i=1:1000)
        i
      
        P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
            polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
        V_MV=P*A_MV_inv;
        V_Be=P*A_Be_inv;

       [k,area_MV] = convhull(V_MV');       
       [k,area_Be] = convhull(V_Be'); 

        samples_poly=double(subs(P*getT(deg,t),t,samples_t))';  
        [k1,area_numeric] = convhull(samples_poly);

        ratios_MV=[ratios_MV area_MV/area_numeric];
        ratios_Be=[ratios_Be area_Be/area_numeric];
    
    end
    
    means_MV=[means_MV mean(ratios_MV)];
    stds_MV=[stds_MV std(ratios_MV)];

    means_Be=[means_Be mean(ratios_Be)];
    stds_Be=[stds_Be std(ratios_Be)];

end

figure; hold on;

color_MV=[123,218,104]/255; %[178,238,166]/255;
color_Be=[179,141,240]/255;

BezierVolumes=shadedErrorBar(all_degs,means_Be,stds_Be,'lineprops',{'-o','Color',color_Be,'markerfacecolor',color_Be,'LineWidth',1},'patchSaturation',0.33);
MinvoVolumes=shadedErrorBar(all_degs,means_MV,stds_MV,'lineprops',{'-o','Color',color_MV,'markerfacecolor',color_MV,'LineWidth',1},'patchSaturation',0.33);
title('\textbf{Area ratio w.r.t. conv($P$), case $k=m=2$}','FontSize',13);xticks(all_degs);

xlabel('$n$');
ylabel('$\frac{\mathrm{area}(\mathrm{conv}(\mathrm{Control\;Points}))}{\mathrm{area}(\mathrm{conv}(P))}$','FontSize',17);
yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); %ylabel('$r$','FontSize', 12);  %'LineSpec','-.'

axes('position',[.2 .4 .5 .4]); box on; tmp=4; xlabel('$n$');
shadedErrorBar(all_degs(1:tmp),means_Be(1:tmp),stds_Be(1:tmp),'lineprops',{'-o','Color',color_Be,'markerfacecolor',color_Be,'LineWidth',1},'patchSaturation',0.33);
yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); %ylabel('$r$','FontSize', 12);  %'LineSpec','-.'
MinvoVolumes=shadedErrorBar(all_degs,means_MV,stds_MV,'lineprops',{'-o','Color',color_MV,'markerfacecolor',color_MV,'LineWidth',1},'patchSaturation',0.33);
axis tight; ylim([0,11]); %xlim([3,4.5])
legend([BezierVolumes.mainLine,MinvoVolumes.mainLine],{'B\''ezier','MINVO'},'Interpreter','latex')

% exportAsPdf(gcf,'area_wrt_convP');

%%%%%%%%%%%%%%%%%%%%Case 3D (k=m=3)

means_MV=[]; means_Be=[];
stds_MV=[]; stds_Be=[];

all_degs=3:7;
for deg=all_degs
    
    ratios_MV=[];
    ratios_Be=[];
    
    A_MV=getA_MV(deg,interv); A_MV_inv=A_MV^(-1);
    A_Be=getA_Be(deg,interv); A_Be_inv=A_Be^(-1);
    tt=linspace(min(interv),max(interv),deg+1);
    for (i=1:1000)
        i
      
        P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
            polyfit(tt,(b-a).*rand(size(tt)) + a,deg);
            polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
        V_MV=P*A_MV_inv;
        V_Be=P*A_Be_inv;

       [k,vol_MV] = convhull(V_MV(1,:),V_MV(2,:),V_MV(3,:));       
       [k,vol_Be] = convhull(V_Be(1,:),V_Be(2,:),V_Be(3,:)); 

        samples_poly=double(subs(P*getT(deg,t),t,samples_t))';  
        [k1,vol_numeric] = convhull(samples_poly(:,1),samples_poly(:,2),samples_poly(:,3));

        ratios_MV=[ratios_MV vol_MV/vol_numeric];
        ratios_Be=[ratios_Be vol_Be/vol_numeric];
    
    end
    means_MV=[means_MV mean(ratios_MV)];
    stds_MV=[stds_MV std(ratios_MV)];
    means_Be=[means_Be mean(ratios_Be)];
    stds_Be=[stds_Be std(ratios_Be)];
end


figure; hold on;

color_MV=[123,218,104]/255; %[178,238,166]/255;
color_Be=[179,141,240]/255;

BezierVolumes=shadedErrorBar(all_degs,means_Be,stds_Be,'lineprops',{'-o','Color',color_Be,'markerfacecolor',color_Be,'LineWidth',1},'patchSaturation',0.33);
MinvoVolumes=shadedErrorBar(all_degs,means_MV,stds_MV,'lineprops',{'-o','Color',color_MV,'markerfacecolor',color_MV,'LineWidth',1},'patchSaturation',0.33);
title('\textbf{Volume ratio w.r.t. conv($P$), case $k=m=3$}','FontSize',13);xticks(all_degs);

xlabel('$n$');
ylabel('$\frac{\mathrm{vol}(\mathrm{conv}(\mathrm{Control\;Points}))}{\mathrm{vol}(\mathrm{conv}(P))}$','FontSize',17);
yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); %ylabel('$r$','FontSize', 12);  %'LineSpec','-.'

axes('position',[.2 .4 .5 .4]); box on; tmp=4; xlabel('$n$');
shadedErrorBar(all_degs(1:tmp),means_Be(1:tmp),stds_Be(1:tmp),'lineprops',{'-o','Color',color_Be,'markerfacecolor',color_Be,'LineWidth',1},'patchSaturation',0.33);
yline(1.0,'-.','Color',[255,159,51]/255,'LineWidth',2); %ylabel('$r$','FontSize', 12);  %'LineSpec','-.'
MinvoVolumes=shadedErrorBar(all_degs,means_MV,stds_MV,'lineprops',{'-o','Color',color_MV,'markerfacecolor',color_MV,'LineWidth',1},'patchSaturation',0.33);
axis tight; ylim([0,11]); %xlim([3,4.5])
legend([BezierVolumes.mainLine,MinvoVolumes.mainLine],{'B\''ezier','MINVO'},'Interpreter','latex')

% exportAsPdf(gcf,'volume_wrt_convP');

