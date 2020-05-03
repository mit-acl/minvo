close all; clear; clc;
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

folder='imgs';
name_figure=[folder,'/minvo_and_bernstein'];

syms t real
interv=[-1,1];
figure;
n_rows=2;
n_cols=4+1;

subplot(n_rows,n_cols,1); hold on
set(gcf, 'Position',  [500, 500, 2000, 550])

T1=[t 1]';
T2=[t*t t 1]';
T3=[t*t*t t*t t 1]';
T4=[t*t*t*t t*t*t t*t t 1]';
T5=[t^5 T4']';
T6=[t^6 T5']';
T7=[t^7 T6']';

for degree=1:4
   
   T=[];
    for i=0:(degree)
       T=[t^i ;T];
    end


   subplot(n_rows,n_cols,degree);
   fplot(getGuessA(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{MINVO, n=',num2str(degree),'}'))
   box on
   
   subplot(n_rows,n_cols,n_cols+degree);
   fplot(computeMatrixForBezier(degree,"m11")*T,interv);
   xlabel('t'); ylim([0,inf]);
   title(strcat('\textbf{Bernstein, n=',num2str(degree),'}'))
   box on
end

% set(gca, 'Position',[0.7813, 0.1100, 0.0371, 0.8150]);

sp_hand1=subplot(n_rows,n_cols,[n_cols, 2*n_cols])
plot(0,0,  0,0,  0,0,  0,0,  0,0,  0,0)
axis off
lgd=legend('$\lambda_1(t)$','$\lambda_2(t)$','$\lambda_3(t)$','$\lambda_4(t)$','$\lambda_5(t)$')
lgd.FontSize = 12;

lgd.Position = [0.9,0.7,1,0.3].*lgd.Position;

pos1 = get(sp_hand1, 'Position') % gives the position of current sub-plot
new_pos1 = [1,1,0.1,1].*pos1 %smaller width
set(sp_hand1, 'Position',new_pos1 ) % set new position of current sub - plot


exportAsPdf(gcf,name_figure)




%%
figure;
subplot(1,2,1);hold on
set(gcf, 'Position',  [500, 500, 3000, 1000])

view1=30;
view2=30;

vx=[0.7 0.5 0 0.1]';
vy=[0.4 1.3 1.1 0.1]';
vz=[1 0 0.8 0]';

A=getGuessA(3,"m11");
pol_x=A'*vx;
pol_y=A'*vy;
pol_z=A'*vz;

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'g');
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

subplot(1,2,2); hold on; 

A=computeMatrixForBezier(3,"m11");
% pol_xbz=A'*vx;
% pol_ybz=A'*vy;
% pol_zbz=A'*vz;

volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A,'b');
fplot3(pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,'r','LineWidth',3);
view(view1, view2)

a1=subplot(1,2,1);
a2=subplot(1,2,2);
allYLim = get([a1 a2], {'YLim'});
allYLim = cat(2, allYLim{:});
set([a1 a2], 'YLim', [min(allYLim), max(allYLim)]);

arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[0 1 0],20,'cylinder',[0.2,0.1]);
arrow3d([0 0 0],[1 0 0],20,'cylinder',[0.2,0.1]);

exportAsSvg(gcf,'imgs/comparison3d')

% text(0,0,1.5,'$x$');

% zlim([min(vz)-0.2,max(vz)+0.2])
% fplot3(pol_x'*T3,pol_y'*T3,0.0*pol_z'*T3,interv,'--r','LineWidth',1);
% fplot3(0.0*pol_x'*T3,pol_y'*T3,pol_z'*T3,interv,4'--r','LineWidth',1);
% fplot3(pol_x'*T3,0.0*pol_y'*T3,pol_z'*T3,interv,'--r','LineWidth',1);

% camlight
% lighting gouraud
% camproj perspective
% view(30, 30)
% shading flat %flat interp
% light
% lighting phong

