close all; clear; clc;

addpath(genpath('./utils'));
addpath(genpath('./solutions'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultfigurecolor',[1 1 1])

interv=[-1,1];

all_deg=2:9;%
all_seg=2:9;

export_figures=false;

for deg=all_deg

figure;

deg

num_col_plots=numel(all_seg);
%%%%
j=1;
for num_seg=all_seg

    %Generate polynomial passing through random points
    a = -1; b = 1;
    tt=linspace(min(interv),max(interv),deg+1);
    P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg); 
        polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];
    
% [vertices_raw, area_union, vertices_union]=myfunction(P, 'MV',interv, num_seg)
subplot(3,num_col_plots,j);
[num_vertices_raw, area_union, area_hull, num_vertices_union]=myfunction(P, 'MV', interv, num_seg);

x_lim_fig=xlim; y_lim_fig=ylim;

%mean(x_lim_fig),y_lim_fig(2)+0.3*(y_lim_fig(2)-y_lim_fig(1))

text(0.5,1.3,['\textbf{num div} \boldmath{$=',num2str(num_seg),'$}'],'HorizontalAlignment','center','Units','normalized')

subplot(3,num_col_plots,num_col_plots+j);
[num_vertices_raw, area_union, area_hull, num_vertices_union]=myfunction(P, 'Be', interv, num_seg);
subplot(3,num_col_plots,2*num_col_plots+j);
[num_vertices_raw, area_union, area_hull, num_vertices_union]=myfunction(P, 'Slefe', interv, num_seg);

j=j+1
end


subplot(3,num_col_plots,1)
x_lim_fig=xlim; y_lim_fig=ylim; ht = text(-0.1,0.5,'\textbf{MINVO}','HorizontalAlignment','center','Units','normalized'); set(ht,'Rotation',90)

subplot(3,num_col_plots,num_col_plots+1)%x_lim_fig(1)-0.1*(x_lim_fig(2)-x_lim_fig(1)),mean(y_lim_fig)
x_lim_fig=xlim; y_lim_fig=ylim; ht = text(-0.1,0.5,'\textbf{B\''ezier}','HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

subplot(3,num_col_plots,2*num_col_plots+1)
x_lim_fig=xlim; y_lim_fig=ylim; ht = text(-0.1,0.5,'\textbf{SLEFE}','HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

set(gcf,'Position',[413         608        1638         400])

if(export_figures)
    exportAsPdf(gcf,['comparison_sleves_deg_',num2str(deg)]);
ebd

end
function [num_vertices_raw, area_union, area_hull, num_vertices_union]=myfunction(P, basis, interv, num_seg)

num_of_breakpts=(num_seg+1);

n=size(P,2)-1; %degree

samples_t=linspace(min(interv),max(interv),num_seg+1);

if(strcmp(basis,'Be'))
    A=getA_Be(n,interv);
    num_vertices_raw=(n+1)*num_seg -  (num_seg-1) ;
    color=[179,141,240]/255;
    name_basis='\textbf{B\''ezier}';

elseif(strcmp(basis,'MV'))
    A=getA_MV(n,interv);
    num_vertices_raw=(n+1)*num_seg;
    color=[123,218,104]/255;
    name_basis='\textbf{MINVO}';
    
elseif(strcmp(basis,'Slefe'))
    
    breakpoints=computeSlefe(P, num_seg, interv);
    color=[255,141,59]/255;
    
    num_vertices_raw=4*num_of_breakpts;
    name_basis='\textbf{Slefe}';
    
else
    error("Not implemented yet")
end

poly_union=[];

for i=1:(length(samples_t)-1)
    if(strcmp(basis,'Be')||strcmp(basis,'MV'))
        a=samples_t(i);
        b=samples_t(i+1);
        P_converted=convertCoeffMatrixFromABtoCD(P,[a,b],interv);
        V=P_converted*inv(A);
    elseif (strcmp(basis,'Slefe'))
        V=[breakpoints{i}.vertices breakpoints{i+1}.vertices];
    end
    [k,av] = convhull(V');

    poly_segment=polyshape(V(1,k),V(2,k)); %V(:,k) are the vertices in the frontier of the convex hull
    
    if(i==1)
        poly_union=poly_segment;
    else
        poly_union=union(poly_union,poly_segment);
    end
    
    
end
 hold on;
plot(poly_union,'FaceColor',color)
area_union=area(poly_union);

%In the lines below, vertices is  [vertex1_x  vertex1_y; 
%                                  vertex2_x  vertex2_y;
%                                  ...]


vertices=filterVerticesPolyShape(poly_union);
plot(vertices(:,1),vertices(:,2),'.k', 'MarkerSize',15)
num_vertices_union=size(vertices,1);
[k,area_hull] = convhull(vertices);

plotCurve(P,interv); %axis off;

% area_union_string=strrep(num2str(area_union,'%.2f'),'+',''); %strrep will change 1e+02--> 1e2 (for shorter notation)
% area_hull_string=strrep(num2str(area_hull,'%.2f'),'+','');   %strrep will change 1e+02--> 1e2 (for shorter notation)

title(['Ve=[',num2str(num_vertices_raw),', ',num2str(num_vertices_union),'] Ar=[',formatNumber(area_union),', ', formatNumber(area_hull),']'],'FontSize',8.2)
box on;
set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca, 'yticklabel',[])

end

function result=formatNumber(number)
    tmp=num2str(number,2);
    if(contains(tmp,'e'))
        %Number is in engineering format (e.g., 2e+5)
        result=strrep(tmp,'+','');
    else
        result=num2str(number,'%.2f');
    end
end


function vertices=filterVerticesPolyShape(poly)

vertices=poly.Vertices;
vertices(any(isnan(vertices), 2), :) = [];
vertices=uniquetol(vertices,'ByRows',true); %See https://www.mathworks.com/help/matlab/ref/uniquetol.html

end

function plotCurve(P,interv)
deg=size(P,2)-1; %degree
syms t real
T=getT(deg,t);  fplot(P(1,:)*T,P(2,:)*T,interv,'r','LineWidth',2);
end

%%


% deg=3;

%Coefficients of the polynomial curve
% P=[ -18.1250    9.3750   31.2500    6.2500  246.8750  274.3750
%     -46.8750  -46.8750  281.2500  -93.7500 -234.3750  180.6250];

% all_vertexes=[];%Its columns are the vertexes
% for i=1:(length(t_breakpoints)-1)
%     a=t_breakpoints(i);
%     b=t_breakpoints(i+1);
%     P_converted=convertCoeffMatrixFromABtoCD(P,[a,b],interv);
%     A_MV=getA_MV(deg,interv); 
%     V=P_converted*inv(A_MV);
%     all_vertexes=[all_vertexes V];
%     %plot_convex_hull(P_converted(1,:)',P_converted(2,:)',P_converted(3,:)',A,'b',0.0017);    
% end
% plot_splitted_convex_hulls(P,A_MV,interv,num_seg,'r',0.01)

% A_MV=getA_MV(deg,interv); 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% A=A_MV;
% num_of_intervals=num_seg;
% color='r';
% radius_sphere=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%
% samples=[];
% 
% all_vertexes=[];%Its columns are the vertexes
% vertices_union=
% plot(poly_union.Vertices(:,1),poly_union.Vertices(:,2),'o')
% figure;  hold on;
%    plot(V(1,k),V(2,k))
%     plot(poly_segment)
%     plot(vertices_segment(1,:),vertices_segment(2,:),'o')
%     
%     all_vertexes=[all_vertexes V];
    %plot_convex_hull(P_converted(1,:)',P_converted(2,:)',P_converted(3,:)',A,'b',0.0017);   
% 
% color_vertex=[.98 .45 .02];
% 
% % for i=1:size(all_vertexes,2)
% %     s1=plotSphere(all_vertexes(:,i),radius_sphere, color_vertex);
% % end
% 
% axis equal
% tmp=gca;
% if (size(findobj(tmp.Children,'Type','Light'))<1) %If still no light in the subplot
%  camlight %create light
% end
% lighting phong
% 
%  
% x=all_vertexes(1,:); y=all_vertexes(2,:); %z=all_vertexes(3,:);
% [k1,area] = convhull(x,y);
% s2=trisurf(k1,x,y,z,'LineWidth',1,'FaceColor',color);
% alpha(s2,0.1)
% 
% k_all=k1(:);
% k_all_unique=unique(k1);
% num_vertexes=length(k_all_unique); %points that are in the frotier of the convex hull
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% % plot_splitted_convex_hulls(P,A,interv,num_of_intervals,'g',0.017);


