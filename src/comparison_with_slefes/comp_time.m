close all; clear; clc;

addpath(genpath('./../utils')); addpath(genpath('./../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])


interv=[-1,1];


all_deg=2:7;%
% all_subdiv=1:6;
all_break_points=2:8;%Only used for Slefe


all_ratios=zeros(numel(all_break_points), numel(all_deg));


for i=1:numel(all_break_points)

    num_seg=all_break_points(i)-1;
    
    i
    
    for j=1:numel(all_deg)
        
        j

        deg=all_deg(j);

        A_MV=getA_MV(deg, interv);
        A_Be=getA_Be(deg, interv);
        A_Be_invA_MV=double(A_Be*inv(A_MV)); %This tabulation is done offline


        tmp_SLEFE=[];
        tmp_MV=[];

        for tmp=1:30

            P=generateRandPol1D(deg,interv); 
            V_Be=P*inv(A_Be);

            [t_break_points, p_down, p_up, comp_time_slefe]=computeSlefeScalarSpeedOptimized(P, deg, num_seg, interv); %Note that, internally, this function uses the Bezier CPs. But the conversion P --> V_Be is not taken into account for the timer
%             [t_break_points2, p_down2, p_up2, comp_time_slefe]=computeSlefeScalar(P, deg, num_seg, interv); 
%         
% %             isequal(p_up,p_up2)
% %             isequal(p_down,p_down2)
            

            tic
            V_Be*A_Be_invA_MV;
            comp_time_MV=toc;

            tmp_SLEFE=[tmp_SLEFE;comp_time_slefe];
            tmp_MV=[tmp_MV;comp_time_MV];
            

        end
        all_ratios(i,j)=mean(tmp_SLEFE)/mean(tmp_MV);
    end
end

vpa(all_ratios,3)


mean(all_ratios(:))