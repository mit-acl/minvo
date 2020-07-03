%This scripts checks numerically the formula of the MINVO paper (volume of
%the convex hull of a polynomial curve)

addpath(genpath('./utils'));

close all; clc; clear;

samples_t=-1:0.01:1;

deg=2;

lateral_areas=[];
J1s=[];
vol_numerics=[];

P=rand(deg,deg+1);

J1=abs(det([P(:,1:end-1)]));

deg_is_odd=(mod(deg,2)==1) ;

if(deg_is_odd)
     a=(deg+1)/2;
     prod=1;
     tt = sym('t', [1,a]);

     for i=0:(a-1)
         for j=0:(a-1)
             if(i>=j)
                 continue;
             end
             ti=tt(i+1); %+1 because of Matlab indexing
             tj=tt(j+1); %+1 because of Matlab indexing
             prod=prod*(ti-tj)^4;
         end
     end

     integrals=prod;
     for i=((deg-1)/2-1):-1:0
          integrals= int( integrals  ,tt(i+2),tt(i+1),1);
     end
     integrals= int( integrals  ,tt(1),-1,1);
     vol_formula=abs(J1*integrals/factorial(deg));

else %deg is even
        
     a=(deg/2)-1;
     tt = sym('t', [1,a+1]);
     
     prod=1;
     
     for i=0:a
         for j=0:a
             if(i>=j)
                 continue;
             end
             ti=tt(i+1); %+1 because of Matlab indexing
             tj=tt(j+1); %+1 because of Matlab indexing
             prod=prod*(ti-tj)^4;
         end
     end
     for i=0:a
         ti=tt(i+1); %+1 because of Matlab indexing
         prod=prod*(ti-1)^2;
     end

     integrals=prod;
     for i=((deg)/2-2):-1:0
         tt(i+1)
         integrals= int( integrals  ,tt(i+2),tt(i+1),1);
     end
     integrals= int( integrals  ,tt(1),-1,1);
%      novale= J1*integrals
     vol_formula=abs(J1*integrals/factorial(deg));
     vpa(vol_formula)
end

syms t
T=[];
 for i=0:deg
     T=[t^i;T];
 end
poly=P*T;        
samples_poly=double(subs(poly,t,samples_t))';  


disp("Computing the volume numerically")
[k1,vol_numeric] = convhulln(samples_poly);
    
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
   J1s=[J1s, J1];
   %lateral_areas=[lateral_areas lateral_area];
    vol_numerics=[vol_numerics vol_numeric];

%%
%     syms alpha beta t0 t1 t2 t3 t4 t5 t6 t real
%     t_vector=(t*ones(deg,1)).^((deg:-1:1)'); %[t^4   t^3  t^2  t];
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
%     ones_vector=ones(deg,1);
    
    
%     if(deg==3)
% %         J2=det([t1^3-t2^3  3*alpha*t1^2  (1-alpha)*3*t2^2;
% %                 t1^2-t2^2    2*alpha*t1  (1-alpha)*2*t2;
% %                   t1-t2        alpha       (1-alpha) ]);
%               
%          J2= det([t1_vector-t2_vector , alpha*t1d_vector  , (1-alpha)*t2d_vector  ]);
%          
%          vol_formula=  abs(vpa(J1*int(int(int(J2,t2, t1,1),t1,-1,1),alpha,0,1)));
%     elseif(deg==5)
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
%     elseif(deg==2)
%         
%         
%         J2=det([ones_vector-t1_vector  ,  (1-alpha)*t1d_vector ] );
%         
% %         J2=det([1-t1^2    (1-alpha)*2*t1 
% %                  1-t1        (1-alpha)   ]);  
%         
%         vol_formula=abs(    J1*int(  int( J2 , t1,-1,1 )      ,alpha,0,1)                   );
%         
%     elseif(deg==4)
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
%      sqrt_deg_g=sqrt(det(g));

% p=patch(samples_poly(:,1),samples_poly(:,2),samples_poly(:,3));
% verts = get(p, 'Vertices');
% faces = get(p, 'Faces');
% a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
% b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));


