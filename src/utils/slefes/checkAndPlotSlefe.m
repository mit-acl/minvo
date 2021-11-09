%This files compares the result produced by our implmentation and the one produced by the "uniexample.c" function of http://www.cise.ufl.edu/research/SurfLab/download/SubLiME.tar.gz

%You simply need to add these lines inside the first loop of the display() function of "uniexample.c":

%     printf("\nVertexs in breakpoint %d\n", i);
%     printf("\nBottomLeft | TopRight \n");
% 
%     printf("%f, %f\n", Lower[i * DIM + 0], Upper[i * DIM + 0]);
%     printf("%f, %f\n", Lower[i * DIM + 1], Upper[i * DIM + 1]);
    
% and this line outside that lopp:

%   printf("\nVertexs in breakpoint %d\n", i);
%   printf("\nBottomLeft | TopRight \n");
% 
%   printf("%f, %f\n", Lower[i * DIM + 0], Upper[i * DIM + 0]);
%   printf("%f, %f\n", Lower[i * DIM + 1], Upper[i * DIM + 1]);
% 
%   printf("\n================\n");

%And then compile and run: `make example` and `./uniexample`

%Note that, once you execute `./uniexample`, you can change the number of
%segments by using the + and - in your keyboard (the max number of segments
%you'll be able to use is defined by MAXSEG (see uniexample.c file)

close all; clear; clc;

addpath(genpath('./../../utils')); addpath(genpath('./../../solutions'));
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1])

V_Be=[
  30.0,  40.0,  0.0;  
  110.0, 340.0, 0.0;  
  250.0, 340.0, 0.0;  
  300.0, 40.0,  0.0;  
  430.0, 40.0,  0.0;  
  550.0, 40.0,  0.0;
]';

deg=5;
num_seg=5;

interv=[-1,1];

P=V_Be*getA_Be(5, interv);

breakpoints=computeSlefe(P, num_seg, interv);

for i=1:numel(breakpoints)
    fprintf("===============\n")
    fprintf("Vertexes in breakpoint %d\n",i-1)
    fprintf("BottomLeft | TopRight\n")
    [sort(unique(double(breakpoints{i}.vertices(1,:))))
     sort(unique(double(breakpoints{i}.vertices(2,:))))]
%     unique(double(slefe{i}.vertices(3,:))) %This is zero (we are in 2D)
end


%% Plot the SLEFE

P=P(1:2,:); %Remove the z coordinate (which is zero)

num_of_breakpts=3;
num_int_subdiv=4;
deg=5;

samples_t=linspace(min(interv),max(interv),num_int_subdiv+1);

figure; hold on;

for i=1:(length(samples_t)-1) %For each one of the subdivisions
    a=samples_t(i);     b=samples_t(i+1); 

    plotCurve(P,[a,b])
    
     breakpoints=computeSlefe(P, num_of_breakpts-1, [a b]);           
     for j=1:(size(breakpoints,2)-1)
         V=double([breakpoints{j}.vertices breakpoints{j+1}.vertices]);
         
         scatter(V(1,:), V(2,:))      
         [k,~] = convhull(V');
         plot(V(1,k),V(2,k))
     end
    
end


function plotCurve(P,interv)
    deg=size(P,2)-1; %degree
    syms t real
    T=getT(deg,t);  fplot(P(1,:)*T,P(2,:)*T,interv,'LineWidth',2);
end