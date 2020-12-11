% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

%This files checks the KKT conditions, see appendix of the paper. 
%These KKT conditions use the Markov–Lukács theorem

clear;  close all; clc;
set(0,'DefaultFigureWindowStyle','docked') %'normal'

addpath(genpath('./utils'));
addpath(genpath('./solutions'));

n=3;
n_is_even = (rem(n, 2) == 0);


syms t real

constraints=[];

G=sym('G_%d%d',[n+1,floor(n/2)+1],'real');
H=sym('H_%d%d',[n+1,floor((n-1)/2)+1],'real');

A=[];
lambdas=[];
for i=1:(n+1)
    i
    coeff_G=G(i,:);
    coeff_H=H(i,:);
    
    TG=[];
    for i=0:(length(coeff_G)-1)
        TG=[t^i;TG]; 
    end

    TH=[];
    for i=0:(length(coeff_H)-1)
        TH=[t^i;TH]; 
    end
    
    g=(coeff_G*TG);
    h=(coeff_H*TH);
    % See now Markov–Lukács theorem   https://www.seas.ucla.edu/~vandenbe/publications/nnp.pdf  Page 952
    if(n_is_even==1) 
        lambdai=(g)^2 + (t+1)*(1-t)*(h)^2; 
    else
        lambdai=(t+1)*(g)^2 + (1-t)*(h)^2;
    end
    coeffs_lambdai=coeffs(lambdai,t,'All');
    A=[A; coeffs_lambdai]; 
    
end


%I want to maximize the absolute value of the determinant of A
detA=computeDet(A);  %I could also write computeDet(A(1:end-1,1:end-1)), because of the constraint A'1=e
% obj=-detA; 
obj=-log(detA*detA); % same as -log(det(A'*A)), and equivalent to minimizing -detA

eq_constraint=sum(A)';
eq_constraint(end)=eq_constraint(end)-1;

l=sym('l',[n+1,1],'real');

lagrangian=obj+l'*eq_constraint;

if(n_is_even)
   error("These KKT conditions below are for odd n")
end


equations=[];
all_sym=[G(:);H(:)];
for i=1:length(all_sym)
    equations=[equations ;diff(lagrangian,all_sym(i))==0 ];
end
equations=[equations ; eq_constraint==0];

% % disp('Solving equations!')
% % s=solve(equations,[C(:);D(:);l(:)]) %Takes forever

n_half=floor(n/2); %=(n-1)/2 when n is odd


ii=2; jj=1;
gi=G(ii,:)';
variable=G(ii,jj);

Rg=[eye(n);zeros(1,n)]+[zeros(1,n);eye(n)];  %Same as toeplitz([1 1 zeros(1,n-1)],[1 zeros(1,n-1)]);
Rh=[-eye(n);zeros(1,n)]+[zeros(1,n);eye(n)]; %Same as toeplitz([-1 1 zeros(1,n-1)],[-1 zeros(1,n-1)]);

znh=zeros(1,n_half);
A_mine=[];
for i=1:(n+1)
    g_i=G(i,:)';
    h_i=H(i,:)';
    AiT=Rg*toeplitz([g_i' znh],[g_i(1) znh])*g_i  + Rh*toeplitz([h_i' znh],[h_i(1) znh])*h_i;
    A_mine=[A_mine;AiT'];  %NOTE THAT i COULD USE conv(a,b) INSTEAD OF toeplitz(f(a))*b (conv also gives multiplication). But in matlab, conv() does not work with symbolic variables. 
end

A_mine2=[rowWiseConv(G,G) rowWiseConv(H,H)]*[Rg';Rh'];

assert(nnz(simplify(A_mine-A))==0);
assert(nnz(simplify(A_mine2-A))==0);

disp("OK, the matrices A coincide!")
%%
LnP1=[zeros(1,n+1);[eye(n) zeros(n,1)]];
Ln=[zeros(1,n);[eye(n-1) zeros(n-1,1)]];

QGij=diff(A,variable);
QGij_mine=2*LnP1^(ii-1)*[  [gi' zeros(1,n-numel(gi))]*((Ln')^(jj-1))*Rg';  
                                        zeros(n,n+1)                 ];

assert(nnz(simplify(QGij-QGij_mine))==0); %Check both results are the same

kkt_eq_mine=trace((-2*inv(A)+l*ones(size(l))')*QGij);

% ee=[zeros(1,n) 1]';
%eq_mine=trace((-2*eye(n+1)+l*ee')*inv(A)*QGij); Note that this also holds
%(because of the constraint A'1=e, but the assert below will not be
%satisfied due to the fact that we are already sustituting

kkt_eq=diff(lagrangian,variable);

assert(nnz(simplify(kkt_eq_mine-kkt_eq))==0);

disp("OK, all tests have passed!")


%%
function result=rowWiseConv(A,B) %MATLAB conv() function doesn't work with symbolic variables

    if(size(A,1)~=size(B,1) || size(A,2)~=size(B,2))
        error("Matrices have to have the same size") 
    end
    
    syms t real
    nn=size(A,2)-1;
    T=getT(nn,t);
    result=[];
    for i=1:size(A,1)
        poly_multiplication=(A(i,:)*T)*(B(i,:)*T); %Note that discrete convolution is the same as polynomial multiplication
        coefficients=coeffs(poly_multiplication,t,'All');
        coefficients=[zeros(1,2*nn-numel(coefficients)+1) coefficients];
        result=[result ;coefficients];
    end
    
end

%% Older code:

% T=toeplitz([gi' znh],[gi(1) znh]);
% should_be=Rg*( diff(T,variable)*gi + T*diff(gi,variable)   )
% mine=2*Rg*(Ln^(jj-1))*[gi;zeros(n-numel(gi),1)]; %Same as Rg*(diff(T,variable)*g1 + T*diff(g1,variable))  
% simplify(should_be-mine)

% syms a b c d f h jj k greal
% syms t real
% Q=sym('Q',[4,4])
% T=getT(3,t);
% tmp=expand(T'*Q*T)
% coeffs(tmp,t,'All')
% 
% 
% tmp=expand(T'*Q*T*t)
% coeffs(tmp,t,'All')

% syms a1 b1 c1 d1 real
% syms a2 b2 c2 d2 real
% syms t real
% f=a1*t^3+b1*t^2+c1*t+d1;
% g=a2*t^3+b2*t^2+c2*t+d2;
% 
% result=f*g;
% coeffs(result,t,'All')'
% 
% To=toeplitz([a1 b1 c1 d1 0 0 0],[a1 0 0 0]);
% To*[a2;b2;c2;d2]
% 
% 
% 
% syms a1 b1 a2 b2 real
% g=a1*t+b1;
% h=a2*t+b2;
% 
% Using t\in [0,1]
% tmp=t*g^2+(1-t)*(h^2)
% should_be=coeffs(tmp,t,'All')'
% result_dos=toeplitz([1 0 0 0],[1 0 0])*toeplitz([a1 b1 0],[a1 0])*[a1;b1]  + toeplitz([-1 1 0 0],[-1 0 0])*toeplitz([a2 b2 0],[a2 0])*[a2;b2];
% 
% simplify(result_uno-result_dos)

% 
% tmp=(t+1)*g^2+(1-t)*(h^2);
% should_be=coeffs(tmp,t,'All')'
% mine_is=toeplitz([1 1 0 0],[1 0 0])*toeplitz([a1 b1 0],[a1 0])*[a1;b1]  + toeplitz([-1 1 0 0],[-1 0 0])*toeplitz([a2 b2 0],[a2 0])*[a2;b2];
% 
% simplify(should_be-mine_is)
% 
% (toeplitz([a1 b1 0],[a1 0])*[a1;b1])'*toeplitz([1 1 0 0],[1 0 0])'
% 
% n=5;
% n_half=floor(n/2);
% g1=sym('g1',[n_half+1,1],'real'); g2=sym('g2',[n_half+1,1],'real'); g3=sym('g3',[n_half+1,1],'real'); g4=sym('g4',[n_half+1,1],'real');g5=sym('g5',[n_half+1,1],'real');g6=sym('g6',[n_half+1,1],'real');
% h1=sym('h1',[n_half+1,1],'real'); h2=sym('h2',[n_half+1,1],'real'); h3=sym('h3',[n_half+1,1],'real'); h4=sym('h4',[n_half+1,1],'real'); h5=sym('h5',[n_half+1,1],'real'); h6=sym('h6',[n_half+1,1],'real');
% 
% 
% tmp=(t+1)*(g1'*getT(numel(g1)-1,t))^2+(1-t)*(h1'*getT(numel(h1)-1,t))^2;
% should_be=coeffs(tmp,t,'All')'
% 
% % m=2;
% nn=2*n_half+1;
% 
% Rg=[eye(nn);zeros(1,n)]+[zeros(1,n);eye(nn)]; %Same as toeplitz([1 1 zeros(1,nn-1)],[1 zeros(1,nn-1)]);
% Rh=[-eye(nn);zeros(1,n)]+[zeros(1,n);eye(nn)]; %Same as toeplitz([-1 1 zeros(1,nn-1)],[-1 zeros(1,nn-1)]);
% 
% mine_is=Rg*toeplitz([g1' zeros(1,n_half)],[g1(1) zeros(1,n_half)])*g1 + Rh*toeplitz([h1' zeros(1,n_half)],[h1(1) zeros(1,n_half)])*h1;
% 
% simplify(should_be-mine_is)
% %%
% znh=zeros(1,n_half);
% 
% A= [g1'*toeplitz([g1' znh],[g1(1) znh])';
%     g2'*toeplitz([g2' znh],[g2(1) znh])';
%     g3'*toeplitz([g3' znh],[g3(1) znh])';
%     g4'*toeplitz([g4' znh],[g4(1) znh])']*Rg'
%    + ...
%    [h1'*toeplitz([h1' znh],[h1(1) znh])';
%     h2'*toeplitz([h2' znh],[h2(1) znh])';
%     h3'*toeplitz([h3' znh],[h3(1) znh])';
%     h4'*toeplitz([h4' znh],[h4(1) znh])']*Rh';
% 
% T=toeplitz([g1' znh],[g1(1) znh]);
% tmp=Rg*T*g1;
% 
% v=g1; ii=1; jj=1; %jj is the element of v
% variable=v(jj);
% 
% 
% should_be=Rg*( diff(T,variable)*v + T*diff(v,variable)   );
% 
% size_L=numel(v)+2;
% L=[zeros(1,size_L);[eye(size_L-1) zeros(size_L-1,1)]];
% mine=2*Rg*(L^(jj-1))*[v;0;0]; %Same as Rg*(diff(T,variable)*g1 + T*diff(g1,variable))  
% 
% simplify(should_be-mine)
% 
% should_be=diff(A,variable);
% mine=2*L^(ii-1)*[[v' 0 0]*(L^(jj-1))'*Rg'; zeros(4,6); ];
% 
% simplify(should_be-mine)

%%
% clc; 
% Q=[A2(1,:)'/norm(A2(1,:)) A2(2,:)'/norm(A2(2,:)) A2(3,:)'/norm(A2(3,:))];
% % http://ee263.stanford.edu/lectures/ellipsoids.pdf    Slide 3
% eigenv1= norm(A2(1,:))^(-2);
% eigenv2= norm(A2(2,:))^(-2);
% eigenv3= norm(A2(3,:))^(-2);
% Lambda=diag([eigenv1 ,eigenv2, eigenv3]);
% 
% S=Q*Lambda*inv(Q)
% 
% tmp=A2'*A2

% tmp=[];
% for i=1:length(all_CDsym)
%     tmp=[tmp ;diff(lagrangian,all_CDsym(i)) ];
% end
% disp('gbasis')
% gbasis(tmp')

% %%
% tmp=A(1:end-1,1:end-1);
% tmp(:,2)=tmp(:,2)-tmp(:,1)
% tmp(:,2)=tmp(:,2)-tmp(:,3)

% 
% 
% n=3;
% n_half=floor(n/2);
% 
% g=sym('g',[n_half+1,1],'real');
% h=sym('h',[n_half+1,1],'real');
% T=getT(n_half,t);
% 
% 
% f=(t+1)*(g'*T)^2  + (1-t)*(h'*T)^2;
% coeff_f=coeffs(f,t,'All')';
% 
% coeff_f_mine_v1=[tSum(g*g'-h*h');0] + [0;tSum(g*g'+h*h')];
% 
% P=[[zeros(1,2*n_half+1) 1] ;[eye(2*n_half+1) zeros(2*n_half+1,1)]];
% I=eye(2*n_half+2);
% g_tilda=[tSum(g*g');0]; h_tilda=[tSum(h*h');0];
% coeff_f_mine_v2=(I+P)*g_tilda-(I-P)*h_tilda;
% 
% assert(double(sum(abs(coeff_f_mine_v1-coeff_f)))<1e-7)
% assert(double(sum(abs(coeff_f_mine_v2-coeff_f)))<1e-7)
% 
% g1=sym('g1',[n_half+1,1],'real'); g2=sym('g2',[n_half+1,1],'real'); g3=sym('g3',[n_half+1,1],'real'); g4=sym('g4',[n_half+1,1],'real');g5=sym('g5',[n_half+1,1],'real');g6=sym('g6',[n_half+1,1],'real');
% h1=sym('h1',[n_half+1,1],'real'); h2=sym('h2',[n_half+1,1],'real'); h3=sym('h3',[n_half+1,1],'real'); h4=sym('h4',[n_half+1,1],'real'); h5=sym('h5',[n_half+1,1],'real'); h6=sym('h6',[n_half+1,1],'real');
% 
% g_tilda1=[tSum(g1*g1');0]; h_tilda1=[tSum(h1*h1');0];
% g_tilda2=[tSum(g2*g2');0]; h_tilda2=[tSum(h2*h2');0];
% g_tilda3=[tSum(g3*g3');0]; h_tilda3=[tSum(h3*h3');0];
% g_tilda4=[tSum(g4*g4');0]; h_tilda4=[tSum(h4*h4');0];
% g_tilda5=[tSum(g5*g5');0]; h_tilda5=[tSum(h5*h5');0];
% g_tilda6=[tSum(g6*g6');0]; h_tilda6=[tSum(h6*h6');0];
% 
% G_hat=[g_tilda1';g_tilda2';g_tilda3';g_tilda4'];
% H_hat=[h_tilda1';h_tilda2';h_tilda3';h_tilda4'];
% 
% A=G_hat*(I+P)'+H_hat*(I-P)'
% 
% % A=[g_tilda1';g_tilda2';g_tilda3';g_tilda4';g_tilda5';g_tilda6']*(I+P)'+[h_tilda1';h_tilda2';h_tilda3';h_tilda4';h_tilda5';h_tilda6']*(I-P)'
% 
% detA=det(A);
% 
% derivate_should_be=diff(detA,g1(1));
% 
% should_be=diff(A,g1(1))
% 
% derivate_mine=det(A)*trace(inv(A)*diff(A,g1(1))) ; %This is good
% 
% simplify(derivate_should_be-derivate_mine)
% 
% %%
% 
% variable=h1(1);
% Q=diff(A,variable);
% 
% % A(2:end,:)=[7 0 0 0; 0 -5 0 0; 0 0 0 20]  (this is to simplify the calculations, if not Matlab fails to simplify correctly sometimes)
% 
% ee=[zeros(n,1);1];
% uno=ones(size(A,1),1);
% l=sym('l',[n+1,1],'real');
% L=log(det(A'*A))+l'*(A'*uno-ee);
% 
% 
% % should_be=diff(log(det(A'*A)),variable);
% % mio=trace(2*inv(A)*Q);
% 
% should_be=diff(L,variable);
% mio=trace((2*inv(A)+l*uno')*Q);
% 
% simplify(expand(should_be-mio))
% 
% % lambda=
% % L=detA*trace(inv(A)*Q)+

% function result=tSum(A)
% 
%     if size(A,1)~=size(A,2)
%         error("Matrix has to be square")
%     end
% 
%     A_flipped=flip(A);
%     result=[];
%     
%     for i=-(size(A,1)-1):1:(size(A,1)-1)
%         result=[result; sum(diag(A_flipped,i))];
%     end
%     
%     
% end
