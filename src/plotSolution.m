% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */
%% Print the ratios of the determinants:
close all; clear; clc;
disp('abs( det(A_MV)/det(A_Be) )')

for i=1:7
    value=abs(det(getSolutionA(i,"m11"))/det(computeMatrixForBezier(i,"m11")));
    vpa(value,4)
end

disp('abs( det(A_MV)/det(A_BS) )')

for i=1:7
    value=abs(det(getSolutionA(i,"m11"))/det(computeMatrixForNonClampedUniformBSpline(i,"m11")));
    vpa(value,4)
end

%%

close all; clear; clc;

addpath(genpath('./utils'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


name_figure='imgs/plots_basis';

syms t real
interv=[-1,1];
figure;
n_rows=8;
n_cols=4;

subplot(n_rows,n_cols,1); hold on
set(gcf, 'Position',  [500, 500, 2000, 2500])

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';

font_size_title=12;

for degree=1:4
   
   T=[];
    for i=0:(degree)
       T=[t^i ;T];
    end


   subplot(n_rows,n_cols,degree);
   fplot(getSolutionA(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{MINVO, n=',num2str(degree),'}'),'FontSize',font_size_title )
   box on
   
   subplot(n_rows,n_cols,n_cols+degree);
   fplot(computeMatrixForBezier(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{Bernstein, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
      
   subplot(n_rows,n_cols,2*n_cols+degree);
   fplot(computeMatrixForNonClampedUniformBSpline(degree,"m11")*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{B-Spline, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,3*n_cols+degree);
   fplot(lagrangePoly(linspace(-1,1,degree+1))*T,interv);
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
   fplot(getSolutionA(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{MINVO, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,5*n_cols+1+(degree-5));
   fplot(computeMatrixForBezier(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{Bernstein, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   
   subplot(n_rows,n_cols,6*n_cols+1+(degree-5));
   fplot(computeMatrixForNonClampedUniformBSpline(degree,"m11")*T,interv);
   xlabel('t'); %ylim([0,inf]);
   title(strcat('\textbf{B-Spline, n=',num2str(degree),'}'),'FontSize',font_size_title)
   box on
   
   subplot(n_rows,n_cols,7*n_cols+1+(degree-5));
   fplot(lagrangePoly(linspace(-1,1,degree+1))*T,interv);
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


exportAsPdf(gcf,name_figure)

%% RESULT for 2D for a given polynomial
figure; hold on;
set(gcf, 'Position',  [500, 500, 3000, 1000])
% subplot(1,2,1);hold on

A=getSolutionA(2,"m11");

v1=[0.1, 0.9];
v2=[0.4  1.0];
v3=[0,  0.35]; 


vx=[v1(1)  v2(1)  v3(1)]';
vy=[v1(2)  v2(2)  v3(2)]';

V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'g','LineWidth',1)
pol_x=A'*vx;
pol_y=A'*vy;


fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);hold on;


%Bernstein
pol_x=pol_x+[0 0 0.5]';
A=computeMatrixForBezier(2,"m11");
vx=inv(A')*pol_x;
vy=inv(A')*pol_y;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'b','LineWidth',1)
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);


%B-Spline
pol_x=pol_x+[0 0 1.5]';
A=computeMatrixForNonClampedUniformBSpline(2,"m11");
vx=inv(A')*pol_x;
vy=inv(A')*pol_y;
v1=[vx(1) vy(1)];
v2=[vx(2) vy(2)];  
v3=[vx(3) vy(3)];
V=[v1; v2; v3];
[k,av] = convhull(V);
fill(V(k,1),V(k,2),'y','LineWidth',1)
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);

alpha 0.2
axis equal
xlim([-2 5.5])
ylim([-1 2.0])

exportAsSvg(gcf,'imgs/comparison2d_poly_given_matlab')

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



A=getSolutionA(3,"m11");
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
A=computeMatrixForBezier(3,"m11");

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2); axis equal;

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);


%%%% BSpline
subplot(1,3,3); hold on; 
A=computeMatrixForNonClampedUniformBSpline(3,"m11");

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

exportAsSvg(gcf,'imgs/comparison3d_poly_given_matlab')

%% RESULT for 2D for a given simplex
figure; hold on;
set(gcf, 'Position',  [500, 500, 3000, 1000])
subplot(1,3,1);hold on

A=getSolutionA(2,"m11");


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

plot(0,0,'*')
fplot(pol_x'*T2,pol_y'*T2,interv,'r','LineWidth',3);
alpha 0.2;xlim([-0.7,0.7]);ylim([-0.1,1.7]);axis equal;

%%%% Bezier
subplot(1,3,2);hold on
A=computeMatrixForBezier(2,"m11");
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
A=computeMatrixForNonClampedUniformBSpline(2,"m11");
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

exportAsSvg(gcf,'imgs/comparison2d_simplex_given_matlab')


%% GEOMETRIC INTERPRETATION OF LAMBDA_I
figure;
subplot(1,2,1);hold on
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


%% RESULT for 3D for a given simplex
figure;
subplot(1,3,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

view1=80;
view2=30;

view1=54;
view2=-4.5;

A=getSolutionA(3,"m11");
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
samples_t=-1:0.01:1;
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
A=computeMatrixForBezier(3,"m11");
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'b',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
axis equal;
ylim([-1.5,1.5]); 

poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=-1:0.01:1;
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
A=computeMatrixForNonClampedUniformBSpline(3,"m11");
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'y',0.02);
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)
axis equal;
ylim([-1.5,1.5]); 

poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
samples_t=-1:0.01:1;
samples_poly=double(subs(poly,t,samples_t));
[k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
camlight
lightangle(gca,45,0)
colormap(winter); axis equal; axis off; 
caxis([0.2 0.7])


% WORKS:
print(gcf,'imgs/comparison3d_simplex_given_matlab','-dpng','-r1000')
 
% DON'T WORK:
% exportAsPdf(gcf,'imgs/comparison_convex_hull')
% saveas(gcf,'imgs/comparison_convex_hull.eps')
% addpath('./utils/plot2svg/plot2svg')
% plot2svg("temperature_standard.svg");
% printeps(get(gcf,'Number'),'imgs/comparison_convex_hull')
% saveas(gcf,'imgs/comparison_convex_hull.png')

%%

% figure; hold on
% volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.015);
% poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
% samples_t=-1:0.01:1;
% samples_poly=double(subs(poly,t,samples_t));
% % centroid_curve=sum(samples_poly,2)/length(samples_t);
% % scatter3(centroid_curve(1),centroid_curve(2),centroid_curve(3),405,'Filled','blue'); 
% % scatter3(samples_poly(1),samples_poly(2),samples_poly(3),405,'Filled','blue'); 
% 
% [k1,av1] = convhull(samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)');
% trisurf(k1,samples_poly(1,:)',samples_poly(2,:)',samples_poly(3,:)','EdgeColor','none','FaceAlpha' ,1.0)%,'FaceColor','cyan'
% fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
% axis equal
% % aplha 0.3
% 
%       camlight
% %      lighting gouraud
%      
%      lightangle(gca,45,0)
% %      lighting gouraud
% %      color_vector=jet;
% %      colormap([color_vector(:,1),color_vector(:,2),color_vector(:,3)])
%      colormap(winter)
%      caxis([0.2 0.7])
%       axis equal
%       axis off
%       
%       view(45, 5)
%       
% % WORKS:
% %  print(gcf,'imgs/comparison_convex_hull_matlab','-dpng','-r1000')
%  
% % DON'T WORK:
% % exportAsPdf(gcf,'imgs/comparison_convex_hull')
% % saveas(gcf,'imgs/comparison_convex_hull.eps')
% % addpath('./utils/plot2svg/plot2svg')
% % plot2svg("temperature_standard.svg");
% % printeps(get(gcf,'Number'),'imgs/comparison_convex_hull')
% % saveas(gcf,'imgs/comparison_convex_hull.png')