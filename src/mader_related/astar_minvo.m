% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


clc; clear; close all; 

q0=[0 0 0]';
q1=[0 0 0]';
q2=[0 0 0]';

Q=[q0 q1 q2];

%Compute upper and lower constraints
%Let's find now v2 (to find q3

qiP1_samples= getChildrenOfQ(Q);

successes=0;

for j=1:size(qiP1_samples,2)
    Q_tmp1=[Q qiP1_samples(:,j)];
    disp("++++++++++++++")
    qiP1_samples1=getChildrenOfQ(Q_tmp1);
    for jj=1:size(qiP1_samples1,2)
        Q_tmp2=[Q_tmp1 qiP1_samples1(:,jj)];
        disp("**")
        qiP1_samples2=getChildrenOfQ(Q_tmp2);
        for jjj=1:size(qiP1_samples1,2)
            Q_tmp3=[Q_tmp2 qiP1_samples1(:,jjj)];
            disp("_")
            qiP1_samples3=getChildrenOfQ(Q_tmp3);
            if(size(qiP1_samples3,2)>0)
                qiP1_samples3;
                successes=successes+1
            end
        end
    end
end

function qiP1_samples=getChildrenOfQ(Q)
    v_max=7;
    p_=3;
    n_samples=3;
    knots_=[0 0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0];
    Mvel_bs2basis=[ 0.5387, 0.08334, -0.03868; 0.5, 0.8333, 0.5;  -0.03867, 0.08333, 0.5387];  
    
    i=size(Q,2);
    
    
    if(mod(i,3)==1)
        Perm=eye(3);
    elseif(mod(i,3)==2)
        Perm=[0 1 0;
              0 0 1;
              1 0 0];
    else
        Perm=[0 0 1;
              1 0 0;
              0 1 0];
    end
        
    Mvel_bs2basis= computeMatrixForBSpline(2,"01")*inv(Perm*getSolutionA(2,"01"));
%     Mvel_bs2basis=Perm*Mvel_bs2basis;
    
%     Mvel_bs2basis=eye(3);
    
    qi=Q(:,end);
    qiM1=Q(:,end-1);
    qiM2=Q(:,end-2);


    viM2=p_ * (qiM1 - qiM2) / (knots_(i - 1 + p_+1) - knots_(i - 1+1)); %one-indexing
    viM1 = p_ * (qi - qiM1) / (knots_(i + p_+1) - knots_(i+1)); %one-indexing 

    Vbs_first_block=[viM2 viM1];

    V_MV=sdpvar(3,3,'full');
    V3_Bs=sdpvar(3,1);                
    constraints=[V_MV==(Vbs_first_block*Mvel_bs2basis(1:2,1:3)+V3_Bs*Mvel_bs2basis(3,1:3))];

    tmp=V_MV(:);
    for j=1:size(tmp,1)
        constraints=[constraints -v_max<=tmp(j)<=v_max];
    end
    obj=sum(V3_Bs);
    result=optimize(constraints,obj,sdpsettings('usex0',0,'solver','fmincon','showprogress',0,'verbose',0,'debug',0 ));
    if(result.problem~=0)
         disp("Problem not solved");
        qiP1_samples=[];
        return;
    end
    V3_Bs_lower= value(V3_Bs);  
    result=optimize(constraints,-obj,sdpsettings('usex0',0,'solver','fmincon','showprogress',0,'verbose',0,'debug',0 ));
    V3_Bs_upper= value(V3_Bs);  

    [Vbs_first_block V3_Bs_upper]*Mvel_bs2basis;
    [Vbs_first_block V3_Bs_lower]*Mvel_bs2basis;


    v2_samples_x=linspace(V3_Bs_lower(1),V3_Bs_upper(1),n_samples);
    v2_samples_y=linspace(V3_Bs_lower(2),V3_Bs_upper(2),n_samples);
    v2_samples_z=linspace(V3_Bs_lower(3),V3_Bs_upper(3),n_samples);
    
%     v2_samples=[linspace(V3_Bs_lower(1),V3_Bs_upper(1),n_samples);
%         linspace(V3_Bs_lower(2),V3_Bs_upper(2),n_samples);
%        linspace(V3_Bs_lower(3),V3_Bs_upper(3),n_samples)];

    qiP1_samples=[];
    for ix=1:size(v2_samples_x,2)
        for iy=1:size(v2_samples_y,2)
            for iz=1:size(v2_samples_z,2)
                vi=[v2_samples_x(ix);v2_samples_y(iy);v2_samples_z(iz)];
                qiP1 = (knots_(i + p_ + 1+1) - knots_(i + 1+1)) * vi / (1.0 * p_) + qi;
                qiP1_samples=[qiP1_samples qiP1];
            end
        end
    end

end


% Qobt=[ -2.5988   -2.5988   -2.5988   -2.5988   -2.5988   -2.5988   -2.5988   -2.5988   -2.5988         0         0         0         0         0         0         0         0         0    2.5988    2.5988    2.5988    2.5988    2.5988    2.5988    2.5988    2.5988;   
%        -2.5988   -2.5988   -2.5988         0         0         0    2.5988    2.5988    2.5988   -2.5988   -2.5988   -2.5988         0         0         0    2.5988    2.5988    2.5988   -2.5988   -2.5988   -2.5988         0         0         0    2.5988    2.5988;
%        -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0    2.5988   -2.5988         0];




