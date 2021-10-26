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

close all; clear; clc;

addpath(genpath('./../utils')); addpath(genpath('./../solutions'));
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

tmp=computeSlefe(P, num_seg, interv);

for i=1:numel(tmp)
    fprintf("===============\n")
    fprintf("Vertexes in breakpoint %d\n",i-1)
    fprintf("BottomLeft | TopRight\n")
    [sort(unique(double(tmp{i}.vertices(1,:))))
     sort(unique(double(tmp{i}.vertices(2,:))))]
%     unique(double(slefe{i}.vertices(3,:))) %This is zero (we are in 2D)
end