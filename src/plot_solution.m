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



%%



set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


name_figure='imgs/plots_basis';

syms t real

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
A=getA_Be(2,interv);
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
A=getA_BS(2,interv);
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

%exportAsSvg(gcf,'imgs/comparison2d_poly_given_matlab')

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

plot(0,0,'*')
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

%exportAsSvg(gcf,'imgs/comparison2d_simplex_given_matlab')


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


%% RESULT for 3D for a given simplex
figure;
subplot(1,3,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

view1=80;
view2=30;

view1=54;
view2=-4.5;

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
%% 
%% RESULT for 3D for a given polynomial


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

function plotAxesArrows(length)
arrow3d([0 0 0],[0 0 length],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 length 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[length 0 0],20,'cylinder',[0.2,0.1]);
end

%print(gcf,'imgs/many_comparisons3d_simplex_given_matlab','-dpng','-r1000')
%%

% figure; hold on
% volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g',0.015);
% poly=[pol_x'*T3,pol_y'*T3,pol_z'*T3]';
% samples_t=min(interv):0.01:max(interv);
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