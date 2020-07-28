clc; clear; close all;
V=rand(2,3);
A=rand(3,4);
A(end,:)= [0 0 0 1]-sum(A(1:end-1,:));

P=V*A;

det_original=det([V' ones(3,1)])

% r= ones(4,1);
r= rand(4,1)%[1 0 0 0]'; %It must be such that det(Aprime!=0)
Aprime=[A ; r'];

ev=[0 0 0 1]';

det_nuevo=(1/det(Aprime))*det([P' ev r]);



A=rand(3,4);
r=rand(4,1);
[A' r]*[V' ones(3,1) zeos(3,1); 0 0 1]