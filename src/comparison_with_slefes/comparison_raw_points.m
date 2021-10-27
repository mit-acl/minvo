%%
close all; clear; clc;

addpath(genpath('./../../utils')); addpath(genpath('./../../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])

figure; 

all_deg=2:7;%
all_subdiv=1:1:6;
all_break_points=all_subdiv+1; %Only used for Slefe

for i=1:numel(all_subdiv)
    
    subplot(1,numel(all_subdiv),i);hold on;
    

    plot(all_deg, getNumRawPointsMV(all_deg,all_subdiv(i)),'-o')
    plot(all_deg, getNumRawPointsBe(all_deg,all_subdiv(i)),'-o')
    
    getNumRawPointsSL(all_deg,all_subdiv(i),3)
    
    for break_points=all_break_points
        h=break_points;
        plot(all_deg, getNumRawPointsSL(all_deg,all_subdiv(i),h),'-o')
    end
    
    
    getNumRawPointsSL(all_deg,all_subdiv(i),3)
    
end


function result=getNumRawPointsMV(n,s)
    result=n*s+s;
end

function result=getNumRawPointsBe(n,s)
    result=n*s+1;
end

function result=getNumRawPointsSL(n,s,h)
    result=4*h*s*ones(size(n));
end