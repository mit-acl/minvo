% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

%This file fits a function to the roots of the MINVO basis functions for n=2,...,7
%The function fitted can be used to estimate the roots of the MINVO basis functions for n>7 
%You may have to run this file several times to get a satisfactory solution. The smallest residual found so far is 0.00501613

close all; clear; clc;
addpath(genpath('./utils'));
addpath(genpath('./solutions'));

interv=[-1,1];

set(0,'DefaultFigureWindowStyle','docked') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(0,'defaultfigurecolor',[1 1 1]);

syms t real

%delete(gcp('nocreate')); %Delete the parallel pool
 
all=[];
 
n_max=7;  figure; hold on;
for n=2:n_max
    subplot(n_max,1,n);  hold on;
    
    rootsMV=getAllRoots_MV(n,interv);
    plot(rootsMV,ones(size(rootsMV)),'o', 'MarkerFaceColor', 'r','Color','r');
    
    inner_roots=rootsMV(rootsMV<1 & rootsMV>-1);
    
%     disp('------------------------------------------------')
    running_index_begin=1;
    running_index_end=1;
    index_cluster=0.0;
    for ii=1:n%numel(idx_unique)
%         disp('-----New cluster')
        running_index_end=running_index_begin+numElCluster(index_cluster,n)-1;
        
        ec=inner_roots(running_index_begin:running_index_end);%    sort(tmp(idx==cluster_index)) %elements in cluster 3
        for j=1:numel(ec) %index inside the cluster
            index_inside_cluster=((j-1));
            index_root=((find(inner_roots==ec(j)))-mean(1:numel(inner_roots)));%/std(1:numel(tmp));
            root_value=ec(j);
            index_basis = getIndexBasisForRoot(root_value, n, interv);
            all=[all; n index_cluster, index_inside_cluster, index_root, index_basis, ec(j)]; %Order is: [degree, index_cluster, index_inside_cluster, index_root, value]
        end
        index_cluster=index_cluster+1;
        running_index_begin=running_index_end+1;
    end
    
    ylabel(['\textbf{n=',num2str(n)','}'])

end

xdata=all(:,1:(end-1));

zdata=all(:,end); %value
b0 = 100*rand(1,3);%-amplitude/2.0;

% [b,resnorm,residual] = lsqcurvefit(@fun,b0,xdata,zdata,[],[],optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt'));
% opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',10000,'StepTolerance',1e-9);
problem = createOptimProblem('lsqcurvefit','x0',b0,'objective',@fun,'lb',[],'ub',[],'xdata',xdata,'ydata',zdata);%,'options',opts
ms = MultiStart('PlotFcns',@gsplotbestf,'Display','iter','UseParallel',true);
[b,errormulti] = run(ms,problem,2550) %b will contain the parameters fitted
%Note that errormulti is norm(fun(b, xdata)-zdata)^2
%%

for n=2:n_max
    figure(1);
    subplot(n_max,1,n);  hold on;
    prediction=fun(b, xdata(xdata(:,1)==n,:));   %[n*ones(n,1),[0:(n-1)]', ]); 
    plot(prediction,ones(size(prediction)),'o', 'Color', 'k', 'MarkerFaceColor', 'k');
end

function result=fun(b,x)

    n=x(:,1);
    ic=x(:,2);%index cluster {0,1,...,(n-1)}
    iic=x(:,3);%index inside cluster {0,1,...,(numElCluster(ic,n)-1)}
    ir=x(:,4);%index root (not used below)
    result=sin((b(1)*(iic-(numElCluster(ic,n)-1)/2.0)    +  b(2)*(ic-(n-1)/2.0))./(n+b(3))); %This is the function that we try to fit

end

function result=numElCluster(ic,n) %ic is the cluster index (0-based), n the degree
   result=floor((n+(isOdd(ic) & isOdd(n)))/2);
end


function result=getIndexBasisForRoot(value_root,n,interv)
    [A rootsA]=getA_MV(n,interv);
    
    for i=1:(n+1)
        if(any(rootsA{i}==value_root))
            result=i-1;
            return
        end
    end
    error("This should never happen")
   
end

% clc
% interv=[-1,1];
% 
% n=26;
% A=getA_MV_Approx(n,interv);
% figure;  fplot(A*getT(n,t),interv)
% 
% %%
% close all;
% figure; fplot(getA_MV(n,interv)*getT(n,t),interv)
% figure; fplot(getA_MV_Approx(n,interv)*getT(n,t),interv)
% 
% getRoots_MV_Approx(7,interv)


% icn=ic-(n-1)/2.0; %index_cluster normalized              This is to enforce symmetry;
% num_elements_cluster=numElCluster(ic,n);%floor((n+(isOdd(ic) & isOdd(n)))/2);
% iicn= x(:,3)-(num_elements_cluster-1)/2.0; %index inside cluster normalized
% centroid=0.0;%sin(1./(b(1)*degree));
% % inside_centroid=iic*b(5) + b(6);
% inside_centroid=sin((b(1)*iic+b(2)*ic)./(degree+b(3)));
% result=centroid+inside_centroid.*abs(1-centroid)^2;

% %%
% hold on;
% 
% for n=2:n_max
%     figure(1);
%     subplot(n_max,1,n);  hold on;
%     prediction=fittedmodel(0:(n-1),n); 
%     plot(prediction,ones(size(prediction)),'o', 'Color', 'k', 'MarkerFaceColor', 'k');
% end

% is_even=~mod(all(:,1),2);
% all=all(~is_even,:)
% n=all_clusters(:,1); %degree
% ic=all_clusters(:,2); %index of the cluster
% z=all_clusters(:,end); %value

% n=2
% dist=1.2*linspace(min(interv),max(interv),n+1)
%  P=lagrangePoly(dist);
%  fplot(P*getT(n,t),interv)
%  
% T1=[t 1]';
% T2=[t*t t 1]';
% T3=[t*t*t t*t t 1]';
% T4=[t*t*t*t t*t*t t*t t 1]';
% T5=[t^5 T4']';
% T6=[t^6 T5']';
% T7=[t^7 T6']';
