%  This file tries to find the smallest triangle that encloses a 2d curve
%  where each coordinate is a third-degree polynomial.
%  Every time you run it, a different random polynomial is used
%  The solver used is Knitro. If it says "Other identified error", run the
%  file again (this happens depending on the initial guess, which is
%  random)


%Tries to find the smallest triangle that encloses the projection of the
%3D MINVO curve (associated with the standard simplex). Note that this
%projection is a rational curve

%Math is in the *.jpg of https://github.com/jtorde/minvo_extra

% rational curve. In this case, we se this polynomial as the p

clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'
addpath(genpath('./../../utils'));
addpath(genpath('./../../solutions'));

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
% constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,2)-1)]   ]; %CHANGED, was sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]
% constraints=[constraints, sum_colums_A(end)==1];


% P=rand(2,4); 
Atmp=getA_MV(3,[-1,1]);
P=[Atmp(2,:);Atmp(4,:)];%For this polynomial...
V=sdpvar(2,3); %I'd like to find the vertexes of the simplex...
obj=-det([V' ones(3,1)],'polynomial');  %that has minimum volume
% obj=det([A(:,1:end-1)],'polynomial')

constraints=[constraints P==V*A];

py=Atmp(3,:)';
constraints=[constraints A'*ones(size(A,1),1)==[0 0 0 1]'-py];

%%%%%%%%%%%%%%%%%%%%%%% INITIAL GUESS
assign(V,rand(2,3))

%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZATION
clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'
% settings=sdpsettings('sparsepop.relaxOrder',3,'savesolveroutput',1,'savesolverinput',1,'solver','sparsepop','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1);
% %,'ipopt.tol',1e-10  %'penlab.max_outer_iter',100000
settings=sdpsettings('usex0',1,'savesolveroutput',0,'savesolverinput',1,'solver','snopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
result=optimize(constraints,obj,settings)

A_minvo=value(A);

syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
    
V=value(V); P=value(P); A=value(A)

figure; hold on
pol=P*T/(1-py'*T);
fplot(pol(1),pol(2), [-1,1])

% plot(subs(pol(1),t,-1/sqrt(3)),subs(pol(2),t,-1/sqrt(3)),'o')
% plot(subs(pol(1),t,1/sqrt(3)),subs(pol(2),t,1/sqrt(3)),'o')

% pol=P*T;
% fplot(pol(1),pol(2), [-1,1])

plot(V(1,:),V(2,:))

[k,av] = convhull(V');
plot(V(1,k),V(2,k))
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1=P(:,1:3);  P2=P(:,4);
% A1=A(:,1:3);  A2=A(:,end);
% A1*inv([P1;-py(1:3)'])*[P2;1-py(4)]
% A2
% 
% q=inv([P1;-py(1:3)'])*[P2;1-py(4)]
% 
% A1'*ones(3,1)+py(1:3)
% figure;
% fplot(A1*[eye(3) q]*T, [-1,1])
% 
% 
% %%
% [V D]=eig(getA_MV(3,[-1,1]))
% for i=1:size(V,2)
% arrow3d([0 0 0],[0 0 1],20,'cylinder',[0.2,0.1]);
% end

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


