%%
close all; clear; clc;

addpath(genpath('./../../utils')); addpath(genpath('./../../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])

figure; 

all_deg=2:7;%
all_subdiv=1:6;
all_break_points=2:8; %Only used for Slefe

    color_Be=[179,141,240]/255;
    color_MV=[123,218,104]/255;
    
    color_SL=flip(colormap(autumn(numel(all_break_points))));

for i=1:numel(all_subdiv)
    
    s=all_subdiv(i);
    handle=subplot(1,numel(all_subdiv)+1,i);hold on; title(['\boldmath{$s=',num2str(s),'$}'])
      
    
    getNumRawPointsSL(all_deg,s,3)
    
    j=1;
    for break_points=all_break_points
        h=break_points;
        plot(all_deg, getNumRawPointsSL(all_deg,s,h),'-o','Color',color_SL(j,:),'MarkerFaceColor',color_SL(j,:),'MarkerEdgeColor','k')
        j=j+1;
    end
   
    plot(all_deg, getNumRawPointsBe(all_deg,s),'-o','Color',color_Be,'MarkerFaceColor',color_Be,'MarkerEdgeColor','k')
    plot(all_deg, getNumRawPointsMV(all_deg,s),'-o','Color',color_MV,'MarkerFaceColor',color_MV,'MarkerEdgeColor','k')
    

    getNumRawPointsSL(all_deg,s,3)
    
    ylim([0,200]); ylabel("\textbf{Num. of raw points}"); xlabel("\boldmath{$n$}"); xticks(all_deg)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% zoom in the image
    if(i==1 || i==2)
        disp("PLOTTING")
        axes('position',get(handle,'position').*[1.0, 1.0, 0.6, 0.25]+[0.02 0.4 0 0]); box on;
        hold on;
        j=1;
        for break_points=all_break_points
            h=break_points;
            plot(all_deg, getNumRawPointsSL(all_deg,s,h),'-o','Color',color_SL(j,:),'MarkerFaceColor',color_SL(j,:),'MarkerEdgeColor','k')
            j=j+1;
        end
        
        plot(all_deg, getNumRawPointsBe(all_deg,s),'-o','Color',color_Be,'MarkerFaceColor',color_Be,'MarkerEdgeColor','k')
        plot(all_deg, getNumRawPointsMV(all_deg,s),'-o','Color',color_MV,'MarkerFaceColor',color_MV,'MarkerEdgeColor','k')
        
        axis tight; box on; xlabel("\boldmath{$n$}"); xticks(all_deg)
        
        if(i==1)
            ylim([0,13.5]);
        else
            ylim([0,26.5]);
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

subplot(1,numel(all_subdiv)+1,numel(all_subdiv)+1);hold on; 

j=1;
for break_points=all_break_points
        h=break_points;
        plot(nan,nan,'-o','Color',color_SL(j,:),'MarkerFaceColor',color_SL(j,:),'MarkerEdgeColor','k')
        j=j+1;
end

plot(nan,nan,'-o','Color',color_Be,'MarkerFaceColor',color_Be,'MarkerEdgeColor','k')
plot(nan, nan,'-o','Color',color_MV,'MarkerFaceColor',color_MV,'MarkerEdgeColor','k')

all_legend={};
for i=1:numel(all_break_points)
all_legend{i}=['\textbf{','SL$_{\mathbf{',num2str(all_break_points(i)),'}}$}'];
end
all_legend{end+1}='\textbf{Be}';
all_legend{end+1}='\textbf{MV}';


lgd=legend(all_legend{:}); 
lgd.Position = [1.01,0.4,1,2.3].*lgd.Position; axis off;

set(gcf,'position',[438         698        1664         282])

exportAsPdf(gcf,'comparison_raw_points');

function result=getNumRawPointsMV(n,s)
    result=n*s+s;
end

function result=getNumRawPointsBe(n,s)
    result=n*s+1;
end

function result=getNumRawPointsSL(n,s,h)
    result=4*h*s*ones(size(n));
end