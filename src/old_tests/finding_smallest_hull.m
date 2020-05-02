clear; clc; close all;
set(0,'DefaultFigureWindowStyle','docked') %'normal'
%Useful to plot the result: http://nurbscalculator.in/
%READ THIS: https://yalmip.github.io/example/nonconvexquadraticprogramming/

W=[];
V=[];
for i=1:4
   W=[W sdpvar(2,2)];
   V=[V sdpvar(2,2)];
end


sdpvar t %This is not optimized, but needed to handle symbolic expressions

constraints=[];

A=[];
sum_Wi=zeros(2,2);
sum_Vi=zeros(2,2);
lambdas=[];
for i=1:4
    Wi=W(:,(2*i-1):2*i);
    Vi=V(:,(2*i-1):2*i);
    sum_Wi=sum_Wi+Wi;
    sum_Vi=sum_Vi+Vi;
    
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

    lambdai=(t+1)*Twi'*Wi*Twi + (1-t)*Tvi'*Vi*Tvi;
    coeffs_lambdai=flip(coefficients(lambdai,t))';
    A=[A; coeffs_lambdai];
        
    
end

sum_colums_A=sum(A);
constraints=[constraints, sum_colums_A(1,1:end-1)==[zeros(1,size(A,1)-1)]   ];
constraints=[constraints, sum_colums_A(end)==1];


a11=A(1,1); a12=A(1,2); a13=A(1,3); a14=A(1,4);
a21=A(2,1); a22=A(2,2); a23=A(2,3); a24=A(2,4);
a31=A(3,1); a32=A(3,2); a33=A(3,3); a34=A(3,4);
a41=A(4,1); a42=A(4,2); a43=A(4,3); a44=A(4,4);

determinant=a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41;

%I want to maximize the absolute value of the determinant of A
obj=-determinant

%I want to maximize the absolute value of the determinant of A
%obj=-abs(det(A,'polynomial')); %Should I put abs() here?

     
A_bezier=[-1 3 -3 1;
          3 -6  3 0;
         -3  3  0 0;
          1  0  0 0]';
   
A_guess=A_bezier;

assign(A,A_guess);

clear t
disp('Starting optimization') %'solver','bmibnb' 'fmincon' ,'solver','sdpt3' 'ipopt' 'knitro' 'scip'
settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'debug',1,'fmincon.maxfunevals',300000,'fmincon.MaxIter', 300000);
% settings=sdpsettings('usex0',1,'savesolveroutput',1,'savesolverinput',1,'solver','fmincon','showprogress',1,'verbose',2,'ipopt.tol',1e-10,'debug',1);
result=optimize(constraints,obj,settings);
%check(constraints)

A_value=value(A);

disp("abs(|A_bezier|/|A_mio|)=")
abs(det(A_bezier)/det(A_value))

%%

 
figure
syms t real
T=[t*t*t t*t t 1]';
fplot(A_value*T, [-1,1])


%%
pol_x=[0.2 0.3 2 1]';%[a b c d]
pol_y=[-0.3 +3 -5 6]';%[a b c d]
pol_z=[1 -0.1 -1 -4]';%[a b c d]



%%
figure; hold on;
fplot3(pol_x'*T,pol_y'*T,pol_z'*T,[-1 1],'r','LineWidth',3);
% axis equal
volumen_mio=plot_convex_hull(pol_x,pol_y,pol_z,A_value,'b');
%volumen_bezier=plot_convex_hull(pol_x,pol_y,pol_z,A_bezier,'g');
disp("abs(|A_bezier|/|A_mio|)=")
abs(det(A_bezier)/det(A_value))
% disp("volumen_mio/volumen_bezier=")
% volumen_mio/volumen_bezier

function constraints=getPsdConstraints(A)
constraints=[];
if(size(A,1)==1)
    constraints=[A(1,1)>=0];
    return;
elseif(size(A,1)==2)
    constraints=[A(1,1)>=0, A(2,2)>=0, (A(1,1)*A(2,2)-A(1,2)*A(2,1))>=0]  
    return;
else
    error('NOT YET IMPLEMENTED')
end
    

    
end

function volume=plot_convex_hull(pol_x,pol_y,pol_z,A,color)
    cx=pol_x;
    cy=pol_y;
    cz=pol_z;

    vx=inv(A')*cx;
    vy=inv(A')*cy;
    vz=inv(A')*cz;

    v1=[vx(1) vy(1) vz(1)]';
    v2=[vx(2) vy(2) vz(2)]';  
    v3=[vx(3) vy(3) vz(3)]';
    v4=[vx(4) vy(4) vz(4)]';  

       
    plot3(v1(1),v1(2),v1(3),'-o','Color',color,'MarkerSize',10)
    plot3(v2(1),v2(2),v2(3),'-o','Color',color,'MarkerSize',10)
    plot3(v3(1),v3(2),v3(3),'-o','Color',color,'MarkerSize',10)
    plot3(v4(1),v4(2),v4(3),'-o','Color',color,'MarkerSize',10)
    
    [k1,volume] = convhull(vx,vy,vz);

    trisurf(k1,vx,vy,vz,'FaceColor',color)
   
    xlabel('x')
    ylabel('y')
    zlabel('z')
    alpha 0.2

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