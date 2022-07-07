%This file compares the volumes/areas obtained by MINVO/Bezier of the curves avaiable in the https://deep-geometry.github.io/abc-dataset
%You need to download a .yml of that ABC dataset first

close all; clc; clear
addpath(genpath('./../utils')); addpath(genpath('./../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])
interv=[-1,1];
syms t;


num_curves_tested=10;
%To be able to use the function ReadYaml, you have to do this:
%Download this: https://code.google.com/archive/p/yamlmatlab/downloads
%Uncompress it and place it in /home/jtorde/installations/YAMLMatlab_0.4.3
%Type edit(fullfile(userpath,'startup.m'))
%Add this line to that file: addpath(genpath('/home/jtorde/installations/YAMLMatlab_0.4.3'))
disp("Reading yaml files")

folder_feat='./abc_dataset/feat/abc_0000_feat_v00/';
folder_obj='./abc_dataset/obj/abc_0000_obj_v00/';

all_files_feat={};
all_files_obj={};

number_parts={'00005045','00004017', '00005022', '00005021', '00002005', '00004021'};

for number_part=number_parts
    tmp=dir([folder_feat,number_part{1},'/']); tmp={tmp.name}; file_feat=[folder_feat,number_part{1},'/',tmp{3}]; %This deletes the . and .. folders
    tmp=dir([folder_obj,number_part{1},'/']); tmp={tmp.name}; file_obj=[folder_obj,number_part{1},'/',tmp{3}]; %This deletes the . and .. folders 
    all_files_feat{end+1}=file_feat;
    all_files_obj{end+1}=file_obj;
end


all_curves=[];
for i=1:numel(all_files_feat)
    file=ReadYaml(all_files_feat{i});
    all_curves=[all_curves, file.curves];
end
disp("Yaml files read")


%%
%Retain only the curves that are BSplines with degree>=3
all_curves_filtered=[];
for i=1:numel(all_curves)
    curve=all_curves{i};
    if(strcmp(curve.type,'BSpline')==1)
        if(curve.degree>=3)
            all_curves_filtered=[all_curves_filtered, curve];
        end
    end
end



% assert(numel(all_curves_filtered)>=num_curves_tested)


close all;

for ii=1:numel(all_curves_filtered)
    nexttile
    if(mod(ii,20)==0)
        figure;
    end
    curve=all_curves_filtered(ii);
      cpoints=cell2mat(curve.poles)';
      knots=cell2mat(curve.knots');      
      spline=spmak(knots,cpoints);
      fnplt(spline,'r',2);
      axis equal
      title(num2str(ii))
end


%Select some of the curves
indexes_filtered=[9,11,18,41,26,61,79,82,106,181];
all_curves_filtered=all_curves_filtered(indexes_filtered);

% for i=indexes_filtered
%     ax=nexttile(i);
%     ax.XColor = [0 0 1];
%     ax.YColor = [0 0 1];
%     ax.ZColor = [0 0 1];
% 
% end

%%
all_deg=3:9;%
all_subdiv=1:6;
all_break_points=2:8;


%AK stands for Average Keeper
all_MV.widthAK(numel(all_subdiv),numel(all_deg))=MyAverageKeeper();
all_Be.widthAK(numel(all_subdiv),numel(all_deg))=MyAverageKeeper();

for j=1:numel(all_break_points)
     all_Slefe{j}.widthAK(numel(all_subdiv),numel(all_deg))=MyAverageKeeper();
end

for ii=1:numel(all_curves_filtered)
%    clc;
   curve=all_curves_filtered(ii);
   curve.type
   if(strcmp(curve.type,'BSpline')==0)
       continue
   end

   if(size( curve.poles,2)~=3)
       error("error")
   end
   
  if(curve.degree>=3)
    
      Q = qr(randn(3)); %Q is a random rotation matrix https://math.stackexchange.com/a/1602779/564801

      disp("Bspline of deg>=3 found")
      cpoints=cell2mat(curve.poles)';
      knots=cell2mat(curve.knots');
      
      knots= (knots - min(knots)) / ( max(knots) - min(knots) ); %normalize between [0,1];
      
      spline=spmak(knots,cpoints);
      sp=fn2fm(spline,'pp'); %Convert to piece-wise polynomial form
      
      %%%
       n_points=100;
       points3D=ppval(sp, linspace(min(sp.breaks), max(sp.breaks), n_points));
%       u=randi([1 numel(sp.breaks)-1],1,1);  %Pick a random segment of the BSpline 
%       points3D=ppval(sp, linspace(sp.breaks(u), sp.breaks(u+1), n_points));
%       for u=randi([1 numel(sp.breaks)-1],1,1)            %1:numel(sp.breaks-1) %For each segment of the Bspline
          %Take sample points
          
          points3D_rotated=Q*points3D;
          points2D=points3D_rotated(1:2,:);

          tt=linspace(min(interv),max(interv),size(points2D,2));
        
          for deg=all_deg
                P= [polyfit(tt,points2D(1,:) ,deg);
                    polyfit(tt,points2D(2,:) ,deg)]
                for num_int=all_subdiv

                    index_num_int=num_int-min(all_subdiv)+1;
                    index_deg=deg-min(all_deg)+1;
                   
                    a=getWidth(P, 'MV', interv, num_int, 0);     all_MV.widthAK(index_num_int,index_deg).addToAverage(a);
                    a=getWidth(P, 'Be', interv, num_int, 0);     all_Be.widthAK(index_num_int,index_deg).addToAverage(a);
                    for j=1:numel(all_break_points)
                         a=getWidth(P, 'Slefe', interv, num_int, all_break_points(j));  all_Slefe{j}.widthAK(index_num_int,index_deg).addToAverage(a);
                    end

                end
          end
    
  end
   
end
%%
for deg=all_deg
    for num_int=all_subdiv
        index_num_int=num_int-min(all_subdiv)+1;
        index_deg=deg-min(all_deg)+1;
        all_MV.width(index_num_int,index_deg)=all_MV.widthAK(index_num_int,index_deg).average;
        all_Be.width(index_num_int,index_deg)=all_Be.widthAK(index_num_int,index_deg).average;

        for j=1:numel(all_break_points)
             all_Slefe{j}.width(index_num_int,index_deg)=all_Slefe{j}.widthAK(index_num_int,index_deg).average;
        end

    end
end

%% Plots
close all;
size_titles=10; nrows=2;  n_col=ceil((size(all_break_points,2)+1)/nrows);

size_numbers_matrix=4.9;
num_plot=1;
%Width
subplot(nrows,n_col,num_plot);        plotMatrix(all_Be.width./all_MV.width, '$r=\frac{[\mathrm{Width}_\mathrm{Be}]}{[\mathrm{Width}_\mathrm{MV}]}$',all_deg, all_subdiv,size_numbers_matrix);
% text(0.5,1.1,'\textbf{Width}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles)
num_plot=num_plot+1;
for j=1:size(all_break_points,2)
    subplot(nrows,n_col,num_plot);  plotMatrix(all_Slefe{j}.width./all_MV.width, ['$r=\frac{[\mathrm{Width}_{\mathrm{SL}_',num2str(all_break_points(j)),'}]}{[\mathrm{Width}_\mathrm{MV}]}$'],all_deg, all_subdiv, size_numbers_matrix);
    num_plot=num_plot+1;
end

num_plot=1;
subplot(nrows,n_col,1); 
ht = text(0.5,1.1,'\textbf{Be vs. MV}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',0.0)
num_plot=num_plot+1;

for j=1:size(all_break_points,2)
    
    subplot(nrows,n_col,num_plot);
    num_plot=num_plot+1;
    ht = text(0.5,1.1,['\textbf{SL$_{\mathbf{',num2str(all_break_points(j)),'}}$ vs. MV}'],'HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',0.0)
    
end
% 
% subplot(nrows,n_col,n_col+1);
% ht = text(-0.3,0.5,'\textbf{SLEFE vs. MINVO}','HorizontalAlignment','center','Units','normalized','FontSize',size_titles); set(ht,'Rotation',90)

set(gcf,'Position',[1032         713        1284         339])
exportAsPdf(gcf,['comparison_slefes_colored_matrices_CAD_width']);

%%
figure;

tiledlayout(2,ceil(numel(all_files_obj)/2));

% j=1;
for file_obj=all_files_obj
     nexttile
%     subplot(1,numel(all_files_obj),j)
    dispObj(readwObj(file_obj{1}));
%     j=j+1;
end


% set(gcf,'Position',[ 649  106  976  1121])
%set(gcf,'Position',[ 812         367        1096         482])
% exportAsPdf(gcf,['comparison_slefes_colored_matrices_CAD_width']);


%%
load('curve_current_reviewer.mat')
all_break_points=2:10;
widths_MV=[];
widths_Be=[];
widths_SL={};
for j=1:numel(all_break_points)
    widths_SL{end+1}=[];
end

for i=1:numel(P)
    widths_MV(end+1)=getWidth(P{i}, 'MV', interv, 1, 0);
    widths_Be(end+1)=getWidth(P{i}, 'Be', interv, 1, 0);
    for j=1:numel(all_break_points)
        widths_SL{j}(end+1)=getWidth(P{i}, 'Slefe', interv, 1, all_break_points(j));
    end
end
clc
vpa(widths_Be./widths_MV,4)
for j=1:numel(all_break_points)
    vpa(widths_SL{j}./widths_MV,4)
end


figure
syms t real;
T8=getT(8,t);
for i=1:numel(P)
    fplot(P{i}(1,:)*T8,P{i}(2,:)*T8,interv,'LineWidth',2); hold on;
end
axis equal;
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
box off; axis off;
exportAsPdf(gcf,['specific_curve_CAD_width_matlab']);

%%

function width=getWidth(P, basis, interv, num_int_subdiv, num_of_breakpts)

    n=size(P,2)-1; %degree
    samples_t=linspace(min(interv),max(interv),num_int_subdiv+1);

    if(strcmp(basis,'Be'))
         A=getA_Be(n,interv); inv_A=inv(A);    
    elseif(strcmp(basis,'MV'))
         A=getA_MV(n,interv); inv_A=inv(A);
    elseif(strcmp(basis,'Slefe'))
            
    else
        error("Not implemented yet")
    end

    width=-inf;
    for i=1:(length(samples_t)-1) %For each one of the subdivisions
        a=samples_t(i);     b=samples_t(i+1); 
        
        if(strcmp(basis,'Be')||strcmp(basis,'MV'))
            % Times:
            %
            %     ├──────────┬──────────────────┬────────────────┤
            % min(interv)    a                  b             max(interv) 
            P_converted=convertCoeffMatrixFromABtoCD(P,[a,b],interv); %This is needed because inv_A in expressed in interv
            V=double(P_converted*inv_A);
            bb=minBoundingBox(V);
            width=max(width,lengthSmallestSideBoundingBox(bb));
            
        elseif (strcmp(basis,'Slefe'))
             breakpoints=computeSlefe(P, num_of_breakpts-1, [a b]);           
             for j=1:(size(breakpoints,2)-1)
                 V=double([breakpoints{j}.vertices breakpoints{j+1}.vertices]);
                 bb=minBoundingBox(V);
                 width=max(width,lengthSmallestSideBoundingBox(bb));
             end
        end               
      
    end
end

%Returns the smallest side of a non-axis-aligned bounding box
function result=lengthSmallestSideBoundingBox(bb)
    result=inf;
    for i=1:size(bb,2)
        for j=1:size(bb,2)
            if(i==j)
                continue
            else
                d=norm(bb(:,i)-bb(:,j));
                if(d<result)
                    result=d;
                end
            end
        end
    end

end


% function area_or_volume=computeAreaOrVolume(points)
% 
% %     points=[x y z]
% 
%     [U,S]=svd(   bsxfun(@minus,points',mean(points')),   0); % See https://www.mathworks.com/matlabcentral/answers/352830-convhull-convhulln-data-is-coplanar-data-is-degenerate-in-at-least-one-dimension-i-know-it-is-but
% 
%     d=diag(S);
%     d=d(d~=0);
% 
%     if(numel(d)==1)
%         error("TODO")
%     elseif(numel(d)==2) %A 2D curve embedded in 3D
%         [~,area_or_volume]=convhull(U*S(:,1:2)); %This is area
%     elseif(numel(d)==3) %A 3D curve embedded in 3D
%         [~,area_or_volume]=convhull(U*S); %This is volume
%     else 
%         error("TODO")
%     end
% 
% end

% 
% [k,area2] = convhull(points(:,1:2));
% [k,area3] = convhull(points);

%   difference= [sumabs(diff(cpoints(:,1)))    sumabs(diff(cpoints(:,2)))    sumabs(diff(cpoints(:,3)))];
%   
%   dimension=nnz(difference); %Indicates if 2D or 3D. 
%   
%   if(dimension==2 && (curve.degree==3 || curve.degree==4))
%       curves2d{end+1}=curve;
%   elseif (dimension==3)
%       curves3d{end+1}=curve;
%   end
%    curves2d{end+1}=