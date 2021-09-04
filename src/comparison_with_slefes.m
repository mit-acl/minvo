close all; clear; clc;

addpath(genpath('./utils'));
addpath(genpath('./solutions'));

set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultfigurecolor',[1 1 1])

interv=[-1,1];

all_deg=2:9;%
all_seg=1:9;

export_figures=false;

for deg=all_deg

figure;

deg

num_col_plots=numel(all_seg);
%%%%
j=1;
for num_seg=all_seg

%Generate polynomial passing through random points
P=generateRandPol(deg,interv);
    
% [vertices_raw, area_union, vertices_union]=myfunction(P, 'MV',interv, num_seg)
subplot(3,num_col_plots,j);
[num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, 'MV', interv, num_seg, true);
x_lim_fig=xlim; y_lim_fig=ylim; text(0.5,1.3,['\textbf{num div} \boldmath{$=',num2str(num_seg),'$}'],'HorizontalAlignment','center','Units','normalized')

subplot(3,num_col_plots,num_col_plots+j);
[num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, 'Be', interv, num_seg, true);
subplot(3,num_col_plots,2*num_col_plots+j);
[num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, 'Slefe', interv, num_seg, true);

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
end

end

%%
% all_MV=[]; all_Be=[]; all_Slefe=[];
% 
% 
% for i=1:1000
%     i
% num_seg=3;
% deg=randi([min(all_deg),max(all_deg)],1);
% P=generateRandPol(deg,interv);
% 
% [a, b, c, d]=myfunction(P, 'MV', interv, num_seg, false); all_MV=[all_MV; [a, b, c, d]];
% [a, b, c, d]=myfunction(P, 'Be', interv, num_seg, false); all_Be=[all_Be; [a, b, c, d]];
% [a, b, c, d]=myfunction(P, 'Slefe', interv, num_seg, false); all_Slefe=[all_Slefe; [a, b, c, d]];
% 
% end

%%
% close all; clc
% figure; hold on;
% x=[all_MV(:,3),all_Be(:,3),all_Slefe(:,3)*0.0];
% boxplot(x,'symbol','','Labels',{'MINVO','B\''ezier','Slefe'}); title('Area hull');
% hAx=gca; hAx.XAxis.TickLabelInterpreter='latex'; hAx.YAxis.TickLabelInterpreter='latex';
% disp('Area hull, Mean, std')
% mean(x)
% std(x)
% 
% figure; hold on;
% x=[all_MV(:,4),all_Be(:,4),all_Slefe(:,4)*0.0];
% boxplot(x,'symbol','','Labels',{'MINVO','B\''ezier','Slefe'}); title('Area union');
% hAx=gca; hAx.XAxis.TickLabelInterpreter='latex'; hAx.YAxis.TickLabelInterpreter='latex';
% disp('Area union, Mean, std')
% mean(x)
% std(x)

%%

all_MV.num_vertices_raw=zeros(numel(all_seg),numel(all_deg));
all_MV.area_union=zeros(numel(all_seg),numel(all_deg));
all_MV.num_vertices_union=zeros(numel(all_seg),numel(all_deg));
all_MV.area_hull=zeros(numel(all_seg),numel(all_deg));
all_MV.num_vertices_hull=zeros(numel(all_seg),numel(all_deg));

all_Be.num_vertices_raw=zeros(numel(all_seg),numel(all_deg));
all_Be.area_union=zeros(numel(all_seg),numel(all_deg));
all_Be.num_vertices_union=zeros(numel(all_seg),numel(all_deg));
all_Be.area_hull=zeros(numel(all_seg),numel(all_deg));
all_Be.num_vertices_hull=zeros(numel(all_seg),numel(all_deg));

counter=1;

for deg=all_deg
    for num_seg=all_seg
        
        counter
        
        tmp_MV=[];    tmp_Be=[];    tmp_Slefe=[];
        for i=1:1
            P=generateRandPol(deg,interv);
            [a, b, c, d, ee]=myfunction(P, 'MV', interv, num_seg, false); tmp_MV=[tmp_MV; [a, b, c, d, ee]];
            [a, b, c, d, ee]=myfunction(P, 'Be', interv, num_seg, false); tmp_Be=[tmp_Be; [a, b, c, d, ee]];
            [a, b, c, d, ee]=myfunction(P, 'Slefe', interv, num_seg, false); tmp_Slefe=[tmp_Slefe; [a, b, c, d, ee]];
        end
        index_num_seg=num_seg-min(all_seg)+1;
        index_deg=deg-min(all_deg)+1;
        
        
        all_MV.num_vertices_raw(index_num_seg,index_deg)=mean(tmp_MV(:,1));
        all_MV.area_union(index_num_seg,index_deg)=mean(tmp_MV(:,2));
        all_MV.num_vertices_union(index_num_seg,index_deg)=mean(tmp_MV(:,3));
        all_MV.area_hull(index_num_seg,index_deg)=mean(tmp_MV(:,4));
        all_MV.num_vertices_hull(index_num_seg,index_deg)=mean(tmp_MV(:,5));
        
        all_Be.num_vertices_raw(index_num_seg,index_deg)=mean(tmp_Be(:,1));
        all_Be.area_union(index_num_seg,index_deg)=mean(tmp_Be(:,2));
        all_Be.num_vertices_union(index_num_seg,index_deg)=mean(tmp_Be(:,3));
        all_Be.area_hull(index_num_seg,index_deg)=mean(tmp_Be(:,4));
        all_Be.num_vertices_hull(index_num_seg,index_deg)=mean(tmp_Be(:,5));
        
        all_Slefe.num_vertices_raw(index_num_seg,index_deg)=mean(tmp_Slefe(:,1));
        all_Slefe.area_union(index_num_seg,index_deg)=mean(tmp_Slefe(:,2));
        all_Slefe.num_vertices_union(index_num_seg,index_deg)=mean(tmp_Slefe(:,3));
        all_Slefe.area_hull(index_num_seg,index_deg)=mean(tmp_Slefe(:,4));
        all_Slefe.num_vertices_hull(index_num_seg,index_deg)=mean(tmp_Slefe(:,5));
        
        counter=counter+1;
    end
end

%%
close all;
n_col=4;
subplot(2,n_col,1);
matrix_value=all_Be.area_union./all_MV.area_union;
imagesc(matrix_value); title ('Area Union'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value); caxis([0 4])

subplot(2,n_col,n_col+1);
matrix_value=all_Slefe.area_union./all_MV.area_union;
imagesc(matrix_value); title ('Area Union'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,2);
matrix_value=all_Be.area_hull./all_MV.area_hull;
imagesc(matrix_value); title ('Area Hull'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,n_col+2);
matrix_value=all_Slefe.area_hull./all_MV.area_hull;
imagesc(matrix_value); title ('Area Hull'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,3);
matrix_value=all_Be.num_vertices_union./all_MV.num_vertices_union;
imagesc(matrix_value); title ('Num Vertices union'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,n_col+3);
matrix_value=all_Slefe.num_vertices_union./all_MV.num_vertices_union;
imagesc(matrix_value); title ('Num Vertices union'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,4);
matrix_value=all_Be.num_vertices_hull./all_MV.num_vertices_hull;
imagesc(matrix_value); title ('Num Vertices hull'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);

subplot(2,n_col,n_col+4);
matrix_value=all_Slefe.num_vertices_hull./all_MV.num_vertices_hull;
imagesc(matrix_value); title ('Num Vertices hull'); doLabelStuff(all_deg, all_seg); plotRectangles(matrix_value);


%%
% subplot(1,4,2)
% imagesc(all_Be.area_union./all_MV.area_union); title ('Area Union');
% doLabelStuff(all_deg, all_seg);
% 
% subplot(1,4,2)
% imagesc(all_Be.num_vertices_raw./all_MV.num_vertices_raw); title ('num_vertices_raw');
% doLabelStuff(all_deg, all_seg);

% for i=1:1000
%     i
% num_seg=3;
% deg=randi([min(all_deg),max(all_deg)],1);
% P=generateRandPol(deg,interv);
% 
% [a, b, c, d]=myfunction(P, 'MV', interv, num_seg, false); all_MV=[all_MV; [a, b, c, d]];
% [a, b, c, d]=myfunction(P, 'Be', interv, num_seg, false); all_Be=[all_Be; [a, b, c, d]];
% [a, b, c, d]=myfunction(P, 'Slefe', interv, num_seg, false); all_Slefe=[all_Slefe; [a, b, c, d]];
% 
% end

%% Functions

function plotRectangles(matrix_value)

hold on;
for i=1:size(matrix_value,1)
    for j=1:size(matrix_value,2)
        if(matrix_value(i,j)>1) %MINVO performs better
            rectangle('Position',[j-0.5,i-0.5,1.0,1.0],'EdgeColor','r')
%             plot([j-0.5,j+0.5],[i-0.5,i+0.5])
            plot(j,i,'.r')
        end
    end
end
end

function doLabelStuff(all_deg, all_seg)
xlabel('\textbf{Degree $n$}'); ylabel('\textbf{Num of intervals $h$}');
c=colorbar; c.Label.String = 'r'; c.Label.Interpreter = 'latex'; c.TickLabelInterpreter= 'latex';
xticklabels(sprintfc('%d',all_deg)); yticklabels(sprintfc('%d',all_seg)); %yticklabels(sprintfc('%d',flip(all_seg)));
set(gca, 'YTick', 1:numel(all_seg)); set(gca, 'XTick', 1:numel(all_deg));
end

function  P=generateRandPol(deg,interv)

a = -1; b = 1;
tt=linspace(min(interv),max(interv),deg+1);
P= [polyfit(tt,(b-a).*rand(size(tt)) + a,deg); 
    polyfit(tt,(b-a).*rand(size(tt)) + a,deg)];

end

function [num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, basis, interv, num_seg, do_plot)

num_of_breakpts=(num_seg+1);

n=size(P,2)-1; %degree

samples_t=linspace(min(interv),max(interv),num_seg+1);

if(strcmp(basis,'Be'))
    A=getA_Be(n,interv); inv_A=inv(A);
    num_vertices_raw=(n+1)*num_seg -  (num_seg-1) ;
    color=[179,141,240]/255;
    name_basis='\textbf{B\''ezier}';

elseif(strcmp(basis,'MV'))
    A=getA_MV(n,interv); inv_A=inv(A);
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
        V=P_converted*inv_A;
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


area_union=area(poly_union);

%In the lines below, vertices is  [vertex1_x  vertex1_y; 
%                                  vertex2_x  vertex2_y;
%                                  ...]


vertices=filterVerticesPolyShape(poly_union);
num_vertices_union=size(vertices,1);
[k,area_hull] = convhull(vertices);
num_vertices_hull=numel(k);

if(do_plot)
    hold on; plot(poly_union,'FaceColor',color)
    plot(vertices(:,1),vertices(:,2),'.k', 'MarkerSize',15)
    plotCurve(P,interv); %axis off;
    title(['Ve=[',num2str(num_vertices_raw),', ',num2str(num_vertices_union),'] Ar=[',formatNumber(area_hull),', ', formatNumber(area_union),']'],'FontSize',8.2)
    box on;
    set(gca,'xtick',[]); set(gca,'xticklabel',[])
    set(gca,'ytick',[]); set(gca, 'yticklabel',[])
end



end

function result=formatNumber(number)
    tmp=num2str(number,2);
    if(contains(tmp,'e'))
        %Number is in engineering format (e.g., 2e+5)
        result=strrep(tmp,'+',''); %strrep will change 1e+02--> 1e2 (for shorter notation)
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


