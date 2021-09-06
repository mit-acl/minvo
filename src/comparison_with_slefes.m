close all; clear; clc;

addpath(genpath('./utils')); addpath(genpath('./solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])


interv=[-1,1];

all_deg=2:9;%
all_seg=1:9;
%% Plot enclosures for random 2D curves
export_figures=true;
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
text(0.5,1.3,[' \boldmath{$h=',num2str(num_seg),'$}'],'HorizontalAlignment','center','Units','normalized')

subplot(3,num_col_plots,num_col_plots+j);
[num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, 'Be', interv, num_seg, true);
subplot(3,num_col_plots,2*num_col_plots+j);
[num_vertices_raw, area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, 'Slefe', interv, num_seg, true);

j=j+1
end


subplot(3,num_col_plots,1)
ht = text(-0.1,0.5,'\textbf{MINVO}','HorizontalAlignment','center','Units','normalized'); set(ht,'Rotation',90)

subplot(3,num_col_plots,num_col_plots+1)%x_lim_fig(1)-0.1*(x_lim_fig(2)-x_lim_fig(1)),mean(y_lim_fig)
ht = text(-0.1,0.5,'\textbf{B\''ezier}','HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

subplot(3,num_col_plots,2*num_col_plots+1)
x_lim_fig=xlim; y_lim_fig=ylim; ht = text(-0.1,0.5,'\textbf{SLEFE}','HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

set(gcf,'Position',[413         608        1638         400])

if(export_figures)
    exportAsPdf(gcf,['comparison_slefes_deg_',num2str(deg)]);
end

end


%%  MonteCarlo analysis for many random polynomialss
%   NOTE: THIS SECTION TAKES ~25 MIN to finish. If you simply want to generate the plots, you can load "data_comparison_slefes.mat" in ./other folder
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
        for i=1:100
            i
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
n_col=4;  size_titles=12;

%Area union
subplot(2,n_col,1);        plotMatrix(all_Be.area_union./all_MV.area_union, '$r=\frac{[\mathrm{Area}_\mathrm{union,\;Be}]}{[\mathrm{Area}_\mathrm{union,\;MV}]}$',all_deg, all_seg);
text(0.5,1.1,'\textbf{Area union}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

subplot(2,n_col,n_col+1);  plotMatrix(all_Slefe.area_union./all_MV.area_union, '$r=\frac{[\mathrm{Area}_\mathrm{union,\;SLEFE}]}{[\mathrm{Area}_\mathrm{union,\;MV}]}$',all_deg, all_seg);

%Area hull
subplot(2,n_col,2);        plotMatrix(all_Be.area_hull./all_MV.area_hull, '$r=\frac{[\mathrm{Area}_\mathrm{hull,\;Be}]}{[\mathrm{Area}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);
text(0.5,1.1,'\textbf{Area hull}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

subplot(2,n_col,n_col+2);  plotMatrix(all_Slefe.area_hull./all_MV.area_hull, '$r=\frac{[\mathrm{Area}_\mathrm{hull,\;SLEFE}]}{[\mathrm{Area}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);

%Num vertices union
subplot(2,n_col,3);        plotMatrix(all_Be.num_vertices_union./all_MV.num_vertices_union, '$r=\frac{[\mathrm{Vert}_\mathrm{union,\;Be}]}{[\mathrm{Vert}_\mathrm{union,\;MV}]}$',all_deg, all_seg);
text(0.5,1.1,'\textbf{Num. vertices union}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

subplot(2,n_col,n_col+3);  plotMatrix(all_Slefe.num_vertices_union./all_MV.num_vertices_union, '$r=\frac{[\mathrm{Vert}_\mathrm{union,\;SLEFE}]}{[\mathrm{Vert}_\mathrm{union,\;MV}]}$',all_deg, all_seg);

%Num vertices hull
subplot(2,n_col,4);        plotMatrix(all_Be.num_vertices_hull./all_MV.num_vertices_hull, '$r=\frac{[\mathrm{Vert}_\mathrm{hull,\;Be}]}{[\mathrm{Vert}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);
text(0.5,1.1,'\textbf{Num. vertices hull}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

subplot(2,n_col,n_col+4);  plotMatrix(all_Slefe.num_vertices_hull./all_MV.num_vertices_hull, '$r=\frac{[\mathrm{Vert}_\mathrm{hull,\;SLEFE}]}{[\mathrm{Vert}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);


subplot(2,n_col,1); 
ht = text(-0.3,0.5,'\textbf{MINVO vs. B\''ezier}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)
subplot(2,n_col,n_col+1);
ht = text(-0.3,0.5,'\textbf{MINVO vs. SLEFE}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)

set(gcf,'Position',[ 475         582        1619         586])

exportAsPdf(gcf,['comparison_slefes_colored_matrices']);

%% Functions

function plotMatrix(matrix_value, label_colorbar, all_deg, all_seg)

imagesc(matrix_value); %title (title_string);

%%%% Label Stuff
xlabel('Degree $n$'); ylabel('Num. of intervals $h$');
c=colorbar;  c.Label.Interpreter = 'latex'; c.Label.String = label_colorbar; c.TickLabelInterpreter= 'latex'; c.Label.FontSize=13;
%set( c.Label,'Rotation',0);; set(c.Label,'HorizontalAlignment','left')

xticklabels(sprintfc('%d',all_deg)); yticklabels(sprintfc('%d',all_seg)); %yticklabels(sprintfc('%d',flip(all_seg)));
set(gca, 'YTick', 1:numel(all_seg)); set(gca, 'XTick', 1:numel(all_deg));
%%%%

%%%Plot rectangles
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
%%%


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

elseif(strcmp(basis,'MV'))
    A=getA_MV(n,interv); inv_A=inv(A);
    num_vertices_raw=(n+1)*num_seg;
    color=[123,218,104]/255;
    
elseif(strcmp(basis,'Slefe'))
    
    breakpoints=computeSlefe(P, num_seg, interv);
    num_vertices_raw=4*num_of_breakpts;
    color=[255,141,59]/255;    
    
else
    error("Not implemented yet")
end

poly_union=[];

for i=1:(length(samples_t)-1)
    if(strcmp(basis,'Be')||strcmp(basis,'MV'))
        a=samples_t(i);     b=samples_t(i+1);
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


vertices_union=filterVerticesPolyShape(poly_union);
num_vertices_union=size(vertices_union,1);
[k,area_hull] = convhull(vertices_union);
num_vertices_hull=numel(unique(k));


assert(num_vertices_hull<=num_vertices_union);
assert(area_hull>=(1-1e-7)*area_union)

if(do_plot)
    hold on; plot(poly_union,'FaceColor',color)
    plot(vertices_union(:,1),vertices_union(:,2),'.k', 'MarkerSize',15)
    plotCurve(P,interv); %axis off;
    title(['Ve=[',num2str(num_vertices_union),', ',num2str(num_vertices_hull),'] Ar=[',formatNumber(area_union),', ', formatNumber(area_hull),']'],'FontSize',8.2)
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
    vertices(any(isnan(vertices), 2), :) = []; %remove the rows with nnan
    vertices=uniquetol(vertices,'ByRows',true); %See https://www.mathworks.com/help/matlab/ref/uniquetol.html

end

function plotCurve(P,interv)
    deg=size(P,2)-1; %degree
    syms t real
    T=getT(deg,t);  fplot(P(1,:)*T,P(2,:)*T,interv,'r','LineWidth',2);
end



