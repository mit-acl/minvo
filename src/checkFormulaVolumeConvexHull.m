% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

%This scripts checks numerically the formula of the MINVO paper (volume of
%the convex hull of a polynomial curve)

addpath(genpath('./utils'));

close all; clc; clear;

samples_t=-1:0.005:1;

n=4;

lateral_areas=[];
J1s=[];
vol_numerics=[];

P=rand(n,n+1);

% J1=abs(det([P(:,1:end-1)]));

n_is_odd=(mod(n,2)==1) ;

disp("***************Checking formula")
prod=1;
for i=0:n
	for j=0:n
        if(i>=j)
            continue;
        end
        prod=prod*(i-j)/(i+j);
    end
end


vol_formula=abs((abs(det([P(:,1:end-1)]))/factorial(n)) * (2^(n*(n+1)/2))*prod); 

syms t
T=[];
 for i=0:n
     T=[t^i;T];
 end
poly=P*T;        
samples_poly=double(subs(poly,t,samples_t))';  
%%

disp("Computing the volume numerically")
if (n>1)
  [k1,vol_numeric] = convhulln(samples_poly);
elseif(n==1)
    vol_numeric=norm( subs(poly,t,1) - subs(poly,t,-1)  );
end
    
if(abs(1-vol_formula/vol_numeric)>1e-3)
    vol_formula;
    vol_numeric;
    vpa(vol_formula/vol_numeric)
    disp("Fomula is not correct");

else
    disp("GOOD")
    %vpa(vol_formula/vol_numeric)
end

   %lateral_area=area3D(samples_poly(:,1),samples_poly(:,2),samples_poly(:,3));
%    J1s=[J1s, J1];
   %lateral_areas=[lateral_areas lateral_area];
    vol_numerics=[vol_numerics vol_numeric];

  
%%
disp("***************Checking formula theorem 15.2 of Geometry of Moment Spaces")
syms t
T=[];
 for i=0:n
     T=[t^i;T];
 end
 
P=[rot90(eye(n)), zeros(n,1)];
P=convertPFrom01toM11(P);


vol_formula=abs((abs(det([P(:,1:end-1)]))/factorial(n)) * (2^(n*(n+1)/2))*prod);

syms k
vol_formula_geometry_moment_spaces=symprod(beta(k,k) , k, 1, n);

vol_formula_geometry_moment_spaces=vpa(vol_formula_geometry_moment_spaces,8);

if(abs(1-vol_formula/vol_formula_geometry_moment_spaces)>1e-3)
    disp("Fomula doesn't match the one from theorem 15.2 of Geometry of Moment Spaces");

else
    disp("GOOD")
    %vpa(vol_formula/vol_numeric)
end


% poly=P*T;        
% samples_poly=double(subs(poly,t,samples_t))';  
% 
% disp("Computing the volume numerically")
% [k1,vol_numeric] = convhulln(samples_poly);
% vol_numeric


function P_converted=convertPFrom01toM11(P)
    
    syms t real
    
    tt=t/2.0+0.5;
    T=[];
    deg=size(P,2)-1;
    for i=0:(deg)
        T=[tt^i T];
    end
    
    tmp=P*T';

    P_converted=[];

    for i=1:size(P,1)
        coefficients=coeffs(tmp(i),'All');
        coefficients=[zeros(1,deg+1-length(coefficients)) coefficients];
        P_converted=[P_converted ;double(vpa(coefficients))];
    end

end

%%

% if(n_is_odd)
%      a=(n+1)/2;
%      prod=1;
%      tt = sym('t', [1,a]);
% 
%      
% % This is without using the closed-form solution of the integrals
% %      for i=0:(a-1)
% %          for j=0:(a-1)
% %              if(i>=j)
% %                  continue;
% %              end
% %              ti=tt(i+1); %+1 because of Matlab indexing
% %              tj=tt(j+1); %+1 because of Matlab indexing
% %              prod=prod*(ti-tj)^4;
% %          end
% %      end
% % 
% %      integrals=prod;
% %      for i=((n-1)/2-1):-1:0
% %           integrals= int( integrals  ,tt(i+2),tt(i+1),1);
% %      end
% %      integrals= int( integrals  ,tt(1),-1,1)
% % %      vpa(integrals)
% %      vol_formula=abs(   (J1/factorial(n))   *   integrals   );
%      
%      
% 
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      prod=1;
%      for i=0:n
%          for j=0:n
%              if(i>=j)
%                  continue;
%              end
%              prod=prod*(i-j)/(i+j);
%          end
%      end
%     
%      vol_formula=(J1/factorial(n)) * (2^(n*(n+1)/2))*prod; 
%      
% 
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      
% 
% else %n is even
%         
%      a=(n/2)-1;
%      tt = sym('t', [1,a+1]);
%      
%      prod=1;
%      
%      for i=0:a
%          for j=0:a
%              if(i>=j)
%                  continue;
%              end
%              ti=tt(i+1); %+1 because of Matlab indexing
%              tj=tt(j+1); %+1 because of Matlab indexing
%              prod=prod*(ti-tj)^4;
%          end
%      end
%      for i=0:a
%          ti=tt(i+1); %+1 because of Matlab indexing
%          prod=prod*(ti-1)^2;
%      end
% 
%      integrals=prod;
%      for i=((n)/2-2):-1:0
%          integrals= int( integrals  ,tt(i+2),tt(i+1),1);
%      end
%      integrals= int( integrals  ,tt(1),-1,1)
% %      vpa(integrals)
% %      novale= J1*integrals
%      vol_formula=abs(J1*integrals/factorial(n));
%      
%      disp("volume should be")
%      vpa(vol_formula)
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%
% 
%      
%      R=zeros(n+2,n+2);
%      R(end,1:end-1)=ones(size(R(end,1:end-1)));
%      R(1:end-1,end)=ones(size(R(1:end-1,end)));
%      
% 
%      
%      prod=1;
%      for i=0:(n)
%          for j=0:(n)       
%              
%              if(i>=j) % || j==(n+1)
%                  continue;
%              end
%              
% %              prod=prod*((ii-jj)^2)/((jj+ii)^2);
%              prod=prod*((i-j))/((j+i));
%          end
%      end
%      
%      
%      vol_closed_form=(J1/factorial(n)) *  2^(n*(n+1)/2)   *prod; 
%      
% %      disp("volume is")
% %      vpa(vol_closed_form)
% %   
%      disp("ratio is")
%      vpa(vol_formula/vol_closed_form)
% %      
%      %%%%%%%%%%%%%%%%%%%%%%%%
%      
%      
% end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%% VERY OLD
%     syms alpha beta t0 t1 t2 t3 t4 t5 t6 t real
%     t_vector=(t*ones(n,1)).^((n:-1:1)'); %[t^4   t^3  t^2  t];
%     td_vector=diff(t_vector,t);
%     
%     t1_vector=subs(t_vector,t,t1);
%     t2_vector=subs(t_vector,t,t2);
%     t3_vector=subs(t_vector,t,t3);
%     t4_vector=subs(t_vector,t,t4);
%     t5_vector=subs(t_vector,t,t5);
%     t6_vector=subs(t_vector,t,t6);
%     
%     t1d_vector=subs(td_vector,t,t1);
%     t2d_vector=subs(td_vector,t,t2);
%     t3d_vector=subs(td_vector,t,t3);
%     t4d_vector=subs(td_vector,t,t4);
%     t5d_vector=subs(td_vector,t,t5);
%     t6d_vector=subs(td_vector,t,t6);    
% 
%         
%     ones_vector=ones(n,1);
    
    
%     if(n==3)
% %         J2=det([t1^3-t2^3  3*alpha*t1^2  (1-alpha)*3*t2^2;
% %                 t1^2-t2^2    2*alpha*t1  (1-alpha)*2*t2;
% %                   t1-t2        alpha       (1-alpha) ]);
%               
%          J2= det([t1_vector-t2_vector , alpha*t1d_vector  , (1-alpha)*t2d_vector  ]);
%          
%          vol_formula=  abs(vpa(J1*int(int(int(J2,t2, t1,1),t1,-1,1),alpha,0,1)));
%     elseif(n==5)
% %         J2=det([t1^5-t3^5   t2^5-t3^5   alpha*5*t1^4    beta*5*t2^4    (1-alpha-beta)*5*t3^4;
% %                 t1^4-t3^4   t2^4-t3^4   alpha*4*t1^3    beta*4*t2^3    (1-alpha-beta)*4*t3^3;
% %                 t1^3-t3^3   t2^3-t3^3   alpha*3*t1^2    beta*3*t2^2     (1-alpha-beta)*3*t3^2;
% %                 t1^2-t3^2   t2^2-t3^2    alpha*2*t1     beta*2*t2       (1-alpha-beta)*2*t3;
% %                   t1-t3      t2-t3         alpha           beta          (1-alpha-beta)  ]);
%          
%          J2=det( [t1_vector-t3_vector  , t2_vector-t3_vector , alpha*t1d_vector, beta*t2d_vector , (1-alpha-beta)*t3d_vector  ]   );
%          vol_formula= abs( J1* int(   int(  int(  int(  int( J2  ,t3,t2,1)     ,t2,t1,1)      ,t1,-1,1)       ,beta,0,1-alpha)       ,alpha,0,1)  );  
%          
% 
%          
%     elseif(n==2)
%         
%         
%         J2=det([ones_vector-t1_vector  ,  (1-alpha)*t1d_vector ] );
%         
% %         J2=det([1-t1^2    (1-alpha)*2*t1 
% %                  1-t1        (1-alpha)   ]);  
%         
%         vol_formula=abs(    J1*int(  int( J2 , t1,-1,1 )      ,alpha,0,1)                   );
%         
%     elseif(n==4)
%                 
%         J2=det([ones_vector-t2_vector ,  t1_vector-t2_vector  ,  beta* t1d_vector   ,   (1-alpha-beta)*t2d_vector     ]);
% %         
% %         J2=det([1-t2^4   t1^4-t2^4   beta*4*t1^3      (1-alpha-beta)*4*t2^3;
% %                 1-t2^3   t1^3-t2^3   beta*3*t1^2      (1-alpha-beta)*3*t2^2;
% %                 1-t2^2   t1^2-t2^2    beta*2*t1        (1-alpha-beta)*2*t2;
% %                 1-t2      t1-t2         beta           (1-alpha-beta)  ]);         
%          
%         vol_formula=abs(   J1*int(  int(  int( int(  J2   ,t2,t1,1)    ,t1,-1,1)    ,beta,0,1-alpha)      ,alpha,0,1)            );
%      
%     else
%         error("Not implemented yet");
%             
%     end

    %%
%figure
%plot(J1s, vol_numerics,'o')

% 
%     J1_m=[P(:,1:end-1)];
%      
%     J2_m=[t1^2-t2^2    2*lambda*t1  (1-lambda)*2*t2;
%             t1-t2        lambda       (1-lambda) ];
%      
%      J=J1_m*J2_m; 
%      
%      g=J'*J;
%      
%      sqrt_n_g=sqrt(det(g));

% p=patch(samples_poly(:,1),samples_poly(:,2),samples_poly(:,3));
% verts = get(p, 'Vertices');
% faces = get(p, 'Faces');
% a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
% b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));



%      %%%%%%%%%%%%%%%%%%%%%%%
%      prod=1;
%      for i=0:(a-1)
%          for j=0:(a-1)
%              if(i>=j)
%                  continue;
%              end
%              ti=tt(i+1); %+1 because of Matlab indexing
%              tj=tt(j+1); %+1 because of Matlab indexing
%              prod=prod*(ti-tj)^4;
%          end
%      end
%      
%      integrals2=prod;
%      for i=((n-1)/2-1):-1:0
%           integrals2= int( integrals2  ,tt(i+2),0,1);
%      end
%      integrals2= int( integrals2  ,tt(1),0,1)
%      
%      integrals2=((2^(a*(2*a-1)))/factorial(a))*integrals2;
%      
%      disp("==================")
%      disp(vpa(integrals2))
%      
%      disp(vpa(integrals))
%      disp("==================")
%      
%      %%%%%%%%%%%%%%%%%%%%%
%      
%      
%      %%%%%%%%%%%%%%%%%%%Another check
%      Khat=zeros(n+1,n+1);
%      for i=0:n
%          for j=0:n
%              if(i==j)
%                  continue
%              end
%              Khat(i+1,j+1)=0.5*(j-i)/(i+j);
%          end
%      end
%      
%      disp("det of Khat is")
%      detKhat=vpa(det(Khat))
% 
% %      prod=1;
% %      for i=0:n
% %          for j=0:n
% %              if(i>=j)
% %                  continue;
% %              end
% %              prod=prod*sqrt(1/2)*(i-j)^(2)/(i+j)^2;
% %          end
% %      end
% 
%      
%      disp("its pf should be")
%      vpa(sqrt(detKhat))
%      %%%%%%%%%%%%%%%%%




     %%%%%%%%%%%%%%%%%%%%%%%CHECK 1
%      prod=1
%      for i=0:a
%          for j=0:a
%              if(i>=j)
%                  continue;
%              end
%              ti=tt(i+1); %+1 because of Matlab indexing
%              tj=tt(j+1); %+1 because of Matlab indexing
%              prod=prod*(ti-tj)^4;
%          end
%      end
%      for i=0:a
%          ti=tt(i+1); %+1 because of Matlab indexing
%          prod=prod*(ti-1)^2;
%      end
% 
%      integrals=prod;
%      for i=((n)/2-2):-1:0
%          tt(i+2)
%          integrals= int( integrals  ,tt(i+2),tt(i+1),1);
%      end
%      integrals= int( integrals  ,tt(1),-1,1);
%      
%      disp("integrals should be")
%      vpa(integrals)
%               
%     
%      integrals=prod;
%      for i=((n)/2-2):-1:0
%          tt(i+2)
%          integrals= int( integrals  ,tt(i+2),0,1);
%      end
%      integrals= int( integrals  ,tt(1),0,1); 
%      
%     integrals= ((2^( 2*a^2+5*a+3   ))/factorial(a+1))*integrals;  
%     
% %      integrals= ((2^(a^2+5*a+3))/factorial(a+1))*integrals;
%      
%      disp("integrals is")
%      vpa(integrals)
%      
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CHECK 1 --> CORRECTO
     
     %%%%%%%%%%%%%%%%%%%Another check
%      Khat=zeros(n+1,n+1);
%      for i=0:n
%          for j=0:n
%              if(i==j)
%                  continue
%              end
%              Khat(i+1,j+1)=(i-j)/(i+j); %0.5*
%          end
%      end
%      
%      disp("det of Khat should be")
%      detKhat=vpa(det(Khat))
%      disp("its pf should be")
%      vpa(sqrt(detKhat))
     


%      disp("det of Khat is")
%      detKhat=vpa(prod)
%      
%      pf=sqrt(prod);
     
     
%      disp("its pf is")
%      vpa(prod/(2^((n+1)/2)))    
%      
     
%      prod=1;
%      for i=0:n
%          for j=0:n
%              if(i>=j)
%                  continue;
%              end
%              prod=prod*(j-i)/(j+i);
%          end
%      end

%       disp("its pf is")
%      vpa(prod/(2^((n+1)/2)))    

     %%%%%%%%%%%%%%%%%
     
%      prod=1;
%      for i=0:n
%          for j=0:n
%              if(i>=j)
%                  continue;
%              end
%              prod=prod*(i-j)/(i+j);
%          end
%      end
     
%      vol_closed_form=(J1/factorial(n)) * ((2^(a^2+5*a+1))/(a+1))*prod; 

%      for i=0:(n)
%          for j=0:(n)
%              if((i+j)==0)
%                  continue;
%              end
%              R(i+1,j+1)=(j-i)/(j+i);
%          end
%      end
%      
%      disp('det_R')
%      det(R)
%  