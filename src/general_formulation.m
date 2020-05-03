clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'
%Useful to plot the result: http://nurbscalculator.in/
%READ THIS: https://yalmip.github.io/example/nonconvexquadraticprogramming/

deg=5;
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
for i=1:(deg+1)
   W=[W sdpvar(size_Wi,size_Wi)];
   V=[V sdpvar(size_Vi,size_Vi)];
end


sdpvar t %This is not optimized, but needed to handle symbolic expressions

constraints=[];

A=[];
lambdas=[];
for i=1:(deg+1)
    Wi=W(:,(i-1)*size_Wi+1:i*size_Wi);
    Vi=V(:,(i-1)*size_Vi+1:i*size_Vi);

    %Wi and Vi are psd matrices <=> All ppal minors are >=0
    constraints=[constraints, getPsdConstraints(Wi), getPsdConstraints(Vi)];   
    
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
constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
constraints=[constraints, sum_colums_A(end)==1];

%I want to maximize the absolute value of the determinant of A
%obj=-computeDet(A); %This is needed when using fmincon (det(A,'polynomial') doesnt work) 
obj=-det(A,'polynomial');

A_bezier=computeMatrixForBezier(deg, "m11");
   
A_guess=A_bezier;

A_guess=getGuessA(deg,"m11");

assign(A,A_guess);

clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'
settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1); %,'ipopt.tol',1e-10
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','ipopt','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
result=optimize(constraints,obj,settings);
%check(constraints)

A_minvo=value(A);

disp("abs(|A_bezier|/|A_minvo|)=")
abs(det(A_bezier)/det(A_minvo))

%%

syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
    

fplot(A_minvo*T, [-1,1])


%%
pol_x=[0.2 0.3 2 1]';%[a b c d]
pol_y=[-0.3 +3 -5 6]';%[a b c d]
pol_z=[1 -0.1 -1 -4]';%[a b c d]



%%
% figure; hold on;
% fplot3(pol_x'*T,pol_y'*T,pol_z'*T,[-1 1],'r','LineWidth',3);
% % axis equal
% volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A_value,'b');
%volumen_bezier=plot_convex_hull(pol_x,pol_y,pol_z,A_bezier,'g');

function constraints=getPsdConstraints(A)
constraints=[];
if(size(A,1)==1)
    constraints=[A(1,1)>=0];
    return;
elseif(size(A,1)==2)
    constraints=[A(1,1)>=0, A(2,2)>=0, (A(1,1)*A(2,2)-A(1,2)*A(2,1))>=0]  
    return;
elseif(size(A,1)==3)

    tmp=A(2:3, 2:3); %delete 1
    constraints=[constraints getPsdConstraints(tmp)] 
    tmp=A([1, 3:end], [1, 3:end]); %delete 2
    constraints=[constraints getPsdConstraints(tmp)]
    tmp=A(1:2, 1:2); %delete 3
    constraints=[constraints getPsdConstraints(tmp)] 
    
    constraints=[constraints computeDet(A)>=0]
    return;
else
    error('NOT YET IMPLEMENTED')
end
    
   
end



% t2=[1 t]';
% lambda1= A_value(:,1)'*T;
% lambda2= A_value(:,2)'*T;
% lambda3= A_value(:,3)'*T;
% lambda4= A_value(:,4)'*T;
% fplot(lambda1,[0,1]); hold on;
% fplot(lambda2,[0,1]);
% fplot(lambda3,[0,1]);
% fplot(lambda4,[0,1]);
% xlim([0 1])
% temporal=t2'*(t*W3 + (1-t)*V3)*t2; %should be lambda1
% coeff_temporal=vpa(coeffs(temporal,t),4)
% 
% disp('wwwwwwwwwwwwwwwwwwwwwwwwww')
% coeff_lambda1=flip(vpa(coeffs(lambda1,t),4)) % [a b c d]
% coeff_lambda2=flip(vpa(coeffs(lambda2,t),4))
% coeff_lambda3=flip(vpa(coeffs(lambda3,t),4))
% coeff_lambda4=flip(vpa(coeffs(lambda4,t),4))


% sum_Vi_value=value(sum_Vi);
% sum_Wi_value=value(sum_Wi);
% C_value=value(C);
% D_value=value(D);
% 
% W1=W_value(:,1:2); V1=V_value(:,1:2);
% W2=W_value(:,3:4); V2=V_value(:,3:4);
% W3=W_value(:,5:6); V3=V_value(:,5:6);
% W4=W_value(:,7:8); V4=V_value(:,7:8);
% 
% t=0.8;
% 
% [1 t]*(t*C_value + D_value)*[1;t]

% px=5*t*t*t+7*t*t+2*t+1;
% py=-4*t*t*t+6*t*t+5*t+6;
% pz=10*t*t*t+2*t*t+3*t+4;

% figure;
% subplot(2,1,1);
% fplot(pol_x'*T,pol_y'*T,[0 1],'r','LineWidth',3)
% xlabel('x'); ylabel('y');
% subplot(2,1,2);
% fplot(pol_x'*T,pol_z'*T,[0 1],'r','LineWidth',3)
% xlabel('x'); zlabel('z');


% coeff_w1=[-32 64 -40 8]/9.0';
% coeff_w2=[64 -112 49 0]/9.0';
% coeff_w3=[-64 80 -17 1]/9.0';
% coeff_w4=[32 -32 8 0]/9.0';
% 
% figure
% fplot(coeff_w1*T,[0,1]); hold on;
% fplot(coeff_w2*T,[0,1]);
% fplot(coeff_w3*T,[0,1]);
% fplot(coeff_w4*T,[0,1]);

% U_value=value(U);
% W_value=value(W);
% V_value=value(V);

% %         vz(4)=v4(3);
%   
%     plot3(v1(1),v1(2),v1(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v2(1),v2(2),v2(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v3(1),v3(2),v3(3),'-o','Color',color,'MarkerSize',10)
%     plot3(v4(1),v4(2),v4(3),'-o','Color',color,'MarkerSize',10)
%   %  end
%     
%          [k1,volume] = convhull(vx,vy,vz);
%  
% %     if color=='b'
%     trisurf(k1,vx,vy,vz,'FaceColor',color)
%    
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     alpha 0.2
% %      end
%     %axis equal
% end


% disp('WWWWWWWWWWWWWWWWWWWWWW')
% constraints_novale=[]
% for i=1:4
%     Wi=W_bezier(:,(2*i-1):2*i);
%     Vi=W_bezier(:,(2*i-1):2*i);
%     
%     Wi(1,1)
%     Vi(1,1)
%     Wi(1,1)*Wi(2,2)-Wi(1,2)*Wi(1,2)
%     Vi(1,1)*Vi(2,2)-Vi(1,2)*Vi(1,2)
%     constraints_novale=[constraints_novale, (Wi(1,1)>=0), (Vi(1,1)>=0)];
%     constraints_novale=[constraints_novale, (Wi(1,1)*Wi(2,2)-Wi(1,2)*Wi(1,2)>=0)];
%     constraints_novale=[constraints_novale, (Vi(1,1)*Vi(2,2)-Vi(1,2)*Vi(1,2)>=0)];
% end
% disp('WWWWWWWWWWWWWWWWWWWWWW')


% W_bezier=[0.0051   -0.0051    3.0000   -3.0003    0.0005   -0.0005    0.0000   -0.0022;
%          -0.0051    0.0051   -3.0003    3.0006   -0.0005    0.0005   -0.0022    1.0043];
% 
% V_bezier=[1.0000   -1.0025    0.0000   -0.0000    0.0000   -0.0003    0.0000   -0.0000
%        -1.0025    1.0051   -0.0000    0.0006   -0.0003    3.0005   -0.0000    0.0043];
% A_guess=10*rand(4,4);

%  A_guess=[ 1  0.5   0.5   1;
%             50  3    0.5   0;
%             3   3       1    -10;
%             10  -30       4   3];
        
% A_guess=ones(4,4);
% W_bezier = ones(2,8);
% V_bezier = ones(2,8);
%      
% W_bezier=[-50   20     4    6     -9    7     15   12;
%           20     1     6    22     7   0.5     12   19]; 
% 
% V_bezier=[0.5     2     -9    2    2    40    20   2;
%            2    100    2   360   40    -5    2   10]; 

% assign(U,A_guess');
% assign(W,W_bezier);
% assign(V,V_bezier);

% A_random=rand(size(A));
% assign(A,A_random);
% assign(U,A_random');
% % W_random=rand(size(W))
% % assign(W,W_random);
% % assign(V,rand(size(V)));

% cx=[5 7 8 6]';
% constraints=constraints(inv(A)*cx==[1 2 3 4]')

%     ui=[Wi(2,2)-Vi(2,2)   ,  -2*Vi(1,2)+Vi(2,2)+2*Wi(1,2) ,  -Vi(1,1)+2*Vi(1,2)+Wi(1,1) , Vi(1,1)]';
%     U=[U; ui'];

% constraints=[constraints A'==U]; %
% C=sum_Wi-sum_Vi;
% D=sum_Vi;



%This constraint are for \sum p_i(t)=1 \forall i
% constraints=[constraints A*ones(4,1)==[0 0 0 1]'];%Sum \lambda_i(t)=1

%These constraints are for \sum p_i(t)=1 \forall i
% constraints=[constraints, C(2,2)==0];
% constraints=[constraints, C(1,2)+C(2,1)+D(2,2)==0];
% constraints=[constraints, C(1,1)+D(1,2)+D(2,1)==0];
% constraints=[constraints, D(1,1)==1];

% constraints=[constraints, A(2,4)==0, A(3,4)==0, A(4,4)==0 , A(3,3)==0 , A(4,3)==0 ];


% A = sdpvar(4,4,'full'); %Should I enforce A symmetric?? NO!! (for Bezier curves, it's symmetric)