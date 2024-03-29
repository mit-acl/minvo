close all; clear; clc;

addpath(genpath('./../utils')); addpath(genpath('./../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])


interv=[-1,1];

all_deg=2:7;%
all_subdiv=1:6;
all_break_points=2:8;%Only used for Slefe %min(all_deg+1):max(all_deg+1); 
smooth=false;

if(smooth==true)
    string_append='_smooth';
else
    string_append='_nonsmooth';
end

rng('shuffle')

%% Plot enclosures for random 2D curves
export_figures=true;
for deg=all_deg

    figure;

    deg

    num_col_plots=numel(all_subdiv);
    %%%%
    j=1;
    for subdiv=all_subdiv
        %Generate polynomial passing through random points
        if(smooth==true)
            P=generateRandSmoothPol(deg,interv);
        else
            P=generateRandPol(deg,interv);
        end


        % [vertices_raw, area_union, vertices_union]=myfunction(P, 'MV',interv, num_int)
        subplot(3,num_col_plots,j);
        [area_union_MV, num_vertices_union_MV, area_hull_MV, num_vertices_hull_MV]=myfunction(P, 'MV', interv, subdiv, 0, true);
        text(0.5,1.3,[' \boldmath{$s=',num2str(subdiv),'$}'],'HorizontalAlignment','center','Units','normalized')

        subplot(3,num_col_plots,num_col_plots+j);
        [area_union_Be, num_vertices_union_Be, area_hull_Be, num_vertices_hull_Be]=myfunction(P, 'Be', interv, subdiv, 0, true);
        subplot(3,num_col_plots,2*num_col_plots+j);
        [area_union_SL, num_vertices_union_SL, area_hull_SL, num_vertices_hull_SL]=myfunction(P, 'Slefe', interv, subdiv, deg+1, true); %Note that we are using deg+1 breakpoints

        if(subdiv>1)
            MV_smaller_union=(area_union_SL>area_union_MV);
            MV_smaller_hull=(area_hull_SL>area_hull_MV);
            if(MV_smaller_union && MV_smaller_hull)
                fprintf("Both, n=%d, s=%d \n", deg, subdiv)
            elseif(MV_smaller_union)
                fprintf("Smaller Union, n=%d, s=%d \n", deg, subdiv)
            elseif(MV_smaller_hull)
                fprintf("Smaller Hull, n=%d, s=%d \n", deg, subdiv)
            end
        end

        j=j+1;
    end


    subplot(3,num_col_plots,1)
    ht = text(-0.1,0.5,'\textbf{MV}','HorizontalAlignment','center','Units','normalized'); set(ht,'Rotation',90)

    subplot(3,num_col_plots,num_col_plots+1)%x_lim_fig(1)-0.1*(x_lim_fig(2)-x_lim_fig(1)),mean(y_lim_fig)
    ht = text(-0.1,0.5,'\textbf{Be}','HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

    subplot(3,num_col_plots,2*num_col_plots+1)
    x_lim_fig=xlim; y_lim_fig=ylim; ht = text(-0.1,0.5,['\textbf{','SL$_{\mathbf{',num2str(deg+1),'}}$}'],'HorizontalAlignment','center','Units','normalized' ); set(ht,'Rotation',90)

    set(gcf,'Position',[413         608        1638         400])

    if(export_figures)
        exportAsPdf(gcf,['comparison_slefes_deg_',num2str(deg),string_append]);
    end

end

disp("Done!")


%%  MonteCarlo analysis for many random polynomials
%   NOTE: THIS SECTION TAKES ~4h to finish. If you simply want to generate the plots, you can load "data_comparison_slefes_smooth.mat" 


all_MV=[];
all_Be=[];
all_Slefe={};


counter=1;

for deg=all_deg
    for num_int=all_subdiv
        
        counter
        
        tmp_MV=[];    tmp_Be=[];    tmp_Slefe=cell(1,size(all_break_points,2));
        for i=1:100
            i
            if(smooth==true)
                P=generateRandSmoothPol(deg,interv);
            else
                P=generateRandPol(deg,interv);
            end
            [a, b, c, d]=myfunction(P, 'MV', interv, num_int, 0, false);     tmp_MV=[tmp_MV;       [a, b, c, d]];
            [a, b, c, d]=myfunction(P, 'Be', interv, num_int, 0, false);     tmp_Be=[tmp_Be;       [a, b, c, d]];
            for j=1:numel(all_break_points)
                [a, b, c, d]=myfunction(P, 'Slefe', interv, num_int, all_break_points(j), false);  tmp_Slefe{j}=[tmp_Slefe{j}; [a, b, c, d]];
            end
            
        end
        index_num_int=num_int-min(all_subdiv)+1;
        index_deg=deg-min(all_deg)+1;
        
        
        all_MV.area_union(index_num_int,index_deg)=            mean(tmp_MV(:,1));
        all_MV.num_vertices_union(index_num_int,index_deg)=    mean(tmp_MV(:,2));
        all_MV.area_hull(index_num_int,index_deg)=             mean(tmp_MV(:,3));
        all_MV.num_vertices_hull(index_num_int,index_deg)=     mean(tmp_MV(:,4));
        
        all_Be.area_union(index_num_int,index_deg)=            mean(tmp_Be(:,1));
        all_Be.num_vertices_union(index_num_int,index_deg)=    mean(tmp_Be(:,2));
        all_Be.area_hull(index_num_int,index_deg)=             mean(tmp_Be(:,3));
        all_Be.num_vertices_hull(index_num_int,index_deg)=     mean(tmp_Be(:,4));
        
        for j=1:numel(all_break_points)
            all_Slefe{j}.area_union(index_num_int,index_deg)=         mean(tmp_Slefe{j}(:,1));
            all_Slefe{j}.num_vertices_union(index_num_int,index_deg)= mean(tmp_Slefe{j}(:,2));
            all_Slefe{j}.area_hull(index_num_int,index_deg)=          mean(tmp_Slefe{j}(:,3));
            all_Slefe{j}.num_vertices_hull(index_num_int,index_deg)=  mean(tmp_Slefe{j}(:,4));
        end
        
        counter=counter+1;
    end
end

%%

%%%%%%%%%%%%%%%%%%%%%
% all_MV.area_union=zeros(numel(all_subdiv),numel(all_deg));
% all_MV.num_vertices_union=zeros(numel(all_subdiv),numel(all_deg));
% all_MV.area_hull=zeros(numel(all_subdiv),numel(all_deg));
% all_MV.num_vertices_hull=zeros(numel(all_subdiv),numel(all_deg));
% 
% all_Be.area_union=zeros(numel(all_subdiv),numel(all_deg));
% all_Be.num_vertices_union=zeros(numel(all_subdiv),numel(all_deg));
% all_Be.area_hull=zeros(numel(all_subdiv),numel(all_deg));
% all_Be.num_vertices_hull=zeros(numel(all_subdiv),numel(all_deg));
% 
% all_Slefe={};
% for j=1:size(all_break_points,2)
%     all_Slefe{j}.area_union=         zeros(numel(all_subdiv),numel(all_deg));
%     all_Slefe{j}.num_vertices_union= zeros(numel(all_subdiv),numel(all_deg));
%     all_Slefe{j}.area_hull=          zeros(numel(all_subdiv),numel(all_deg));
%     all_Slefe{j}.num_vertices_hull=  zeros(numel(all_subdiv),numel(all_deg));
% end
%%%%%%%%%%%%%%%%%%%%%



close all;
n_col=4;  size_titles=10; nrows=size(all_break_points,2)+1;
size_numbers_matrix=4.6;
%Area union
subplot(nrows,n_col,1);        plotMatrix(all_Be.area_union./all_MV.area_union, '$r=\frac{[\mathrm{Area}_\mathrm{union,\;Be}]}{[\mathrm{Area}_\mathrm{union,\;MV}]}$',all_deg, all_subdiv, size_numbers_matrix);
text(0.5,1.1,'\textbf{Area union}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

for j=1:size(all_break_points,2)
    subplot(nrows,n_col,j*n_col+1);  plotMatrix(all_Slefe{j}.area_union./all_MV.area_union, ['$r=\frac{[\mathrm{Area}_{\mathrm{union,\;SL}_',num2str(all_break_points(j)),'}]}{[\mathrm{Area}_\mathrm{union,\;MV}]}$'],all_deg, all_subdiv, size_numbers_matrix);
end


%Area hull
subplot(nrows,n_col,2);        plotMatrix(all_Be.area_hull./all_MV.area_hull, '$r=\frac{[\mathrm{Area}_\mathrm{hull,\;Be}]}{[\mathrm{Area}_\mathrm{hull,\;MV}]}$',all_deg, all_subdiv, size_numbers_matrix);
text(0.5,1.1,'\textbf{Area hull}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

% subplot(nrows,n_col,n_col+2);  plotMatrix(all_Slefe.area_hull./all_MV.area_hull, '$r=\frac{[\mathrm{Area}_\mathrm{hull,\;SLEFE}]}{[\mathrm{Area}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);

for j=1:size(all_break_points,2)
    subplot(nrows,n_col,j*n_col+2);  plotMatrix(all_Slefe{j}.area_hull./all_MV.area_hull, ['$r=\frac{[\mathrm{Area}_{\mathrm{hull,\;SL}_',num2str(all_break_points(j)),'}]}{[\mathrm{Area}_\mathrm{hull,\;MV}]}$'],all_deg, all_subdiv, size_numbers_matrix);
end


%Num vertices union
subplot(nrows,n_col,3);        plotMatrix(all_Be.num_vertices_union./all_MV.num_vertices_union, '$r=\frac{[\mathrm{Vert}_\mathrm{union,\;Be}]}{[\mathrm{Vert}_\mathrm{union,\;MV}]}$',all_deg, all_subdiv, size_numbers_matrix);
text(0.5,1.1,'\textbf{Num. vertices union}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

% subplot(nrows,n_col,n_col+3);  plotMatrix(all_Slefe.num_vertices_union./all_MV.num_vertices_union, '$r=\frac{[\mathrm{Vert}_\mathrm{union,\;SLEFE}]}{[\mathrm{Vert}_\mathrm{union,\;MV}]}$',all_deg, all_seg);

for j=1:size(all_break_points,2)
    subplot(nrows,n_col,j*n_col+3);  plotMatrix(all_Slefe{j}.num_vertices_union./all_MV.num_vertices_union, ['$r=\frac{[\mathrm{Vert}_{\mathrm{union,\;SL}_',num2str(all_break_points(j)),'}]}{[\mathrm{Vert}_\mathrm{union,\;MV}]}$'],all_deg, all_subdiv, size_numbers_matrix);
end

%Num vertices hull
subplot(nrows,n_col,4);        plotMatrix(all_Be.num_vertices_hull./all_MV.num_vertices_hull, '$r=\frac{[\mathrm{Vert}_\mathrm{hull,\;Be}]}{[\mathrm{Vert}_\mathrm{hull,\;MV}]}$',all_deg, all_subdiv, size_numbers_matrix);
text(0.5,1.1,'\textbf{Num. vertices hull}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)

% subplot(nrows,n_col,n_col+4);  plotMatrix(all_Slefe.num_vertices_hull./all_MV.num_vertices_hull, '$r=\frac{[\mathrm{Vert}_\mathrm{hull,\;SLEFE}]}{[\mathrm{Vert}_\mathrm{hull,\;MV}]}$',all_deg, all_seg);

for j=1:size(all_break_points,2)
    subplot(nrows,n_col,j*n_col+4);  plotMatrix(all_Slefe{j}.num_vertices_hull./all_MV.num_vertices_hull, ['$r=\frac{[\mathrm{Vert}_{\mathrm{hull,\;SL}_',num2str(all_break_points(j)),'}]}{[\mathrm{Vert}_\mathrm{hull,\;MV}]}$'],all_deg, all_subdiv, size_numbers_matrix);
end

subplot(nrows,n_col,1); 
ht = text(-0.3,0.5,'\textbf{Be vs. MV}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)


for j=1:size(all_break_points,2)
    
    subplot(nrows,n_col,j*n_col+1);
    
    ht = text(-0.3,0.5,['\textbf{SL$_{\mathbf{',num2str(all_break_points(j)),'}}$ vs. MV}'],'HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)
    
end
% 
% subplot(nrows,n_col,n_col+1);
% ht = text(-0.3,0.5,'\textbf{SLEFE vs. MINVO}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)

% set(gcf,'Position',[ 649  106  976  1121])
set(gcf,'Position',[ 649  106  976  1170])
exportAsPdf(gcf,['comparison_slefes_colored_matrices',string_append]);

%% Functions
function [area_union, num_vertices_union, area_hull, num_vertices_hull]=myfunction(P, basis, interv, num_int_subdiv, num_of_breakpts, do_plot)

n=size(P,2)-1; %degree

samples_t=linspace(min(interv),max(interv),num_int_subdiv+1);

if(strcmp(basis,'Be'))
    A=getA_Be(n,interv); inv_A=inv(A);
    color=[179,141,240]/255;

elseif(strcmp(basis,'MV'))
    A=getA_MV(n,interv); inv_A=inv(A);
    color=[123,218,104]/255;
    
elseif(strcmp(basis,'Slefe'))
    color=[255,141,59]/255;    
else
    error("Not implemented yet")
end



poly_union=[];

for i=1:(length(samples_t)-1) %For each one of the subdivisions
    a=samples_t(i);     b=samples_t(i+1); 
    
    if(strcmp(basis,'Be')||strcmp(basis,'MV'))
        % Times:
        %
        %     ├──────────┬──────────────────┬────────────────┤
        % min(interv)    a                  b             max(interv) 
        P_converted=convertCoeffMatrixFromABtoCD(P,[a,b],interv); %This is needed because inv_A in expressed in interv
        V=double(P_converted*inv_A);
        [k,~] = convhull(V');
        poly_segment=polyshape(V(1,k),V(2,k)); %V(:,k) are the vertices in the frontier of the convex hull
        if(isempty(poly_union))  poly_union=poly_segment; else poly_union=union(poly_union,poly_segment); end
        
    elseif (strcmp(basis,'Slefe'))
         breakpoints=computeSlefe(P, num_of_breakpts-1, [a b]);           
         for j=1:(size(breakpoints,2)-1)
             V=double([breakpoints{j}.vertices breakpoints{j+1}.vertices]);
             [k,~] = convhull(V');
             poly_segment=polyshape(V(1,k),V(2,k)); %V(:,k) are the vertices in the frontier of the convex hull
             if(isempty(poly_union))  poly_union=poly_segment; else poly_union=union(poly_union,poly_segment); end
         end
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
    title(['Ve=[',num2str(num_vertices_union),', ',num2str(num_vertices_hull),'] Ar=[',formatNumber(area_union,3),', ', formatNumber(area_hull,3),']'],'FontSize',8.2)
    box on;
    set(gca,'xtick',[]); set(gca,'xticklabel',[])
    set(gca,'ytick',[]); set(gca, 'yticklabel',[])
end



end







