% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'

deg=3;

tmp=sdpvar(deg, deg+1);
A=[tmp; [zeros(1,deg) 1]-sum(tmp)]; %With this I can avod adding the constraints sum_colums_A=[0 0 0 1];
% A=sdpvar(deg+1, deg+1, 'full');

constraints=[];
% sum_colums_A=sum(A);
% constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
% constraints=[constraints, sum_colums_A(end)==1];

A_solution=getSolutionA(deg,"m11"); 

syms tt real
T=[];
for i=0:(deg)
   T=[tt^i ;T];
end

roots_lambda_i=getRootsLambdai(deg,"m11"); 
% 
%%%%%%%%%%%%%%%% Impose the positivity constraints on lambda_i on the roots
%of lambda_i
% for i=1:size(A_solution,1)
%        
%     tmp=roots_lambda_i{i}
%     
%     for j=1:length(tmp)
%       T_j=double(subs(T,tt,tmp(j))); 
%       constraints=[constraints, A(i,:)*T_j>=0]; %Impose this constraint for each of the roots of lambda_i
%     end      
%     
% end
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Impose the positivity constraints on lambda_i on ALL the
% roots of the lambdas
all_roots=[];
for i=1:size(A_solution,1)
    all_roots=[all_roots roots_lambda_i{i}];
end

for i=1:size(A_solution,1)
    for j=1:length(all_roots)
      T_j=double(subs(T,tt,all_roots(j))); 
      constraints=[constraints, A(i,:)*T_j>=0]; %Impose this constraint for all the roots
    end      
    
end

constraints


%%%%%%%%%%%%%%%
A_cropped=A(1:end-1,1:end-1);
obj=-det(A_cropped,'polynomial'); % abs(det(A))==abs(det(A_cropped)) because of the constraint sum_rows_A=[0 0 ... 0 1]

assign(A,getSolutionA(deg,"m11"))

% check(constraints)

clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip' %,'ipopt.tol',1e-10  %'penlab.max_outer_iter',100000

%%

settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% settings=sdpsettings('moment.order',4,'savesolveroutput',1,'savesolverinput',1,'solver','moment','showprogress',1,'verbose',2,'debug',1);
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1);
% settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
result=optimize(constraints,obj,settings)

relaxation_order=3;
if(deg<=3)
    relaxation_order=3; %WORKS
end

% [sol,x,momentdata] = solvemoment(constraints,obj,[],relaxation_order);


% check(constraints)

A_minvo=value(A);

disp("abs(|A_minvo|)=")
abs(det(A_minvo))

%%

syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
    

fplot(A_minvo*T, [-1,1])
