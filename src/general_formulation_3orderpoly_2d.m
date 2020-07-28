%  This file tries to find the smallest triangle that encloses a 2d curve
%  where each coordinate is a third-degree polynomial.
%  Every time you run it, a different random polynomial is used
%  The solver used is Knitro. If it says "Other identified error", run the
%  file again (this happens depending on the initial guess, which is
%  random)

% This file allowed me to check that the optimal matrix A for this case
% changes for each polynomial. 

clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'

deg=3; %right now only works for deg=3
deg_is_even = (rem(deg, 2) == 0);

if(deg_is_even==1)
    d=deg/2;
    size_Wi=d+1;
    size_Vi=d;
else
    d=(deg-1)/2;
    size_Wi=d+1;
    size_Vi=d+1;
end


W=[]; V=[];
for i=1:(deg)  %CHANGED, was i=1:(deg+1)
   W=[W sdpvar(size_Wi,size_Wi)];
   V=[V sdpvar(size_Vi,size_Vi)];
end

%%
sdpvar t %This is not optimized, but needed to handle symbolic expressions

constraints=[];

A=[];
lambdas=[];
for i=1:(deg) %CHANGED, was i=1:(deg+1)
    Wi=W(:,(i-1)*size_Wi+1:i*size_Wi);
    Vi=V(:,(i-1)*size_Vi+1:i*size_Vi);

    %Wi and Vi are psd matrices <=> All ppal minors are >=0
    constraints=[constraints, getPsdConstraints(Wi), getPsdConstraints(Vi)];   
%     constraints=[constraints, Wi>=0, Vi>=0]; %Works when using solvemoment sdp
    
    Tvi=[];
    for i=0:(size(Vi,1)-1)
        Tvi=[Tvi; t^i]; 
    end

    Twi=[];
    for i=0:(size(Wi,1)-1)
        Twi=[Twi; t^i]; 
    end

    if(deg_is_even==1)
        lambdai=Twi'*Wi*Twi + (t+1)*(1-t)*Tvi'*Vi*Tvi;
    else
        lambdai=(t+1)*Twi'*Wi*Twi + (1-t)*Tvi'*Vi*Tvi;
    end
    coeffs_lambdai=flip(coefficients(lambdai,t))';
    A=[A; coeffs_lambdai]; 
    
end


sum_colums_A=sum(A); 
constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,2)-1)]   ]; %CHANGED, was sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]
constraints=[constraints, sum_colums_A(end)==1];


P=rand(2,4); %For this polynomial...
V=sdpvar(2,3); %I'd like to find the vertexes of the simplex...
obj=-det([V' ones(3,1)],'polynomial');  %that has minimum volume

constraints=[constraints P==V*A];

%%%%%%%%%%%%%%%%%%%%%%% INITIAL GUESS
assign(V,rand(2,3))

%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION
clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'
% settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1);
% %,'ipopt.tol',1e-10  %'penlab.max_outer_iter',100000
settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','knitro','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
result=optimize(constraints,obj,settings)

A_minvo=value(A);

syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
    
V=value(V); P=value(P); A=value(A)
fplot(A_minvo*T, [-1,1])
figure; hold on
pol=P*T;
fplot(pol(1),pol(2), [-1,1])
plot(V(1,:),V(2,:))

[k,av] = convhull(V');
plot(V(1,k),V(2,k))


function constraints=getPsdConstraints(A)
constraints=[];

%See Theorem 17 Of Lecture 5 of https://homepages.cwi.nl/~monique/eidma-seminar-parrilo/Parrilo-LectureNotes-EIDMA.pdf

sdpvar lambda;
n=size(A,1);
charac_poly=computeDet(lambda*eye(size(A,1))-A);

coeffs_charac_poly=coefficients(charac_poly,lambda); %[p0 p1 ... p_n-1 1]

size(coeffs_charac_poly)
for(i=0:(n-1))
    j=i+1;
    constraints=[constraints coeffs_charac_poly(j)*((-1)^(n-i))>=0 ];
end

   
end


