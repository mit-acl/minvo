close all; clear; clc;

addpath(genpath('./utils')); addpath(genpath('./solutions'));

figure; hold on;  syms t real;

deg=2;
n_interv=5;
T=getT(deg,t);

all_V=[ 1 0.5 1 3 3.25 4.25 5 
        -1 0.25 3 0 3.25 1 2 ]

% colors=['k','r','c','g','m'];
default_colors = colororder(); 

n_int_knots=n_interv-1;
deltaT=1/(n_int_knots+1);
interm=deltaT*(1:n_int_knots);
knots = [zeros(1,deg+1)   interm   (max(interm)+deltaT)*ones(1,deg+1)];

interv=[min(knots), min(knots)+deltaT]

for i=1:n_interv
    segment_key=i-1;

    A=computeMatrixForAnyBSpline(deg, deg+1+segment_key, knots, interv)

    V=all_V(:,i:i+deg);
    
    P=V*A;
    pol_x=P(1,:)'; pol_y=P(2,:)';
    subplot(1,2,1); hold on;
    patch(V(1,:),V(2,:),default_colors(i,:),'FaceAlpha',.1)
    fplot(pol_x'*T,pol_y'*T,interv,'LineWidth',3, 'Color',default_colors(i,:));
    plot(V(1,:), V(2,:), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k')

    i

    subplot(1,2,2); hold on;
    
    tmp_x=[min(interv), max(interv), max(interv), min(interv)];
    tmp_y=[0, 0, 1, 1];

    patch(tmp_x,tmp_y,default_colors(i,:),'FaceAlpha',.1)

    fplot(A(1,:)*T,interv,'k','LineWidth',1);
    fplot(A(2,:)*T,interv,'k','LineWidth',1);
    fplot(A(3,:)*T,interv,'k','LineWidth',1);
    ylim([0,1])
    xline(max(interv))
    xlabel('$t$')

    interv=interv+deltaT;

end




subplot(1,2,1); hold on;

xlabel('$x(t)$')
ylabel('$y(t)$')

subplot(1,2,2);
title('Basis functions')
xlabel('$t$')

set(gcf, 'Position',  [603         589        1890         300])
exportAsPdf(gcf,'spline_example')



%     if(segment_key<= n_interv-(deg-1))
%         disp("First part")
%         A=computeMatrixForClampedUniformBSpline(deg, segment_key, interv);
%     else
%         disp("Second part")
%         A=computeMatrixForClampedUniformBSpline(deg, segment_key-n_interv+1, interv);
%     end


% 
% V=[ 0.1000    1.0000    0.8000  ;
%     0.5000    0.2000    1.0000  ];


