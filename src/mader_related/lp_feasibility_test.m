% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */


clc; clear; close all;  

    %For MINVO getA_MV
    %For Bezier  getA_Be
    
    Mvel_bs2basis0= computeMatrixForClampedUniformBSpline(2,0,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis1= computeMatrixForClampedUniformBSpline(2,1,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis2= computeMatrixForClampedUniformBSpline(2,2,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis3= computeMatrixForClampedUniformBSpline(2,3,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis4= computeMatrixForClampedUniformBSpline(2,4,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis5= computeMatrixForClampedUniformBSpline(2,5,"01")*inv(getA_MV(2,"01"))
    Mvel_bs2basis6= computeMatrixForClampedUniformBSpline(2,5,"01")*inv(getA_MV(2,"01"))

    v_max=7.0;
    viM2=[0 0 0]';
    viM1=[8.5 0 0]';
        
    constraints=[];
    %%%%%%%%%%%%%%% INTERVAL 0
    a=sdpvar(3,1);  
    b=sdpvar(3,1); 
    c=sdpvar(3,1);
    d=sdpvar(3,1);
    e=sdpvar(3,1);
    
    tmp=[viM2 viM1 a]*Mvel_bs2basis1;
    tmp=tmp(:);
    for j=1:length(tmp)
        constraints=[constraints -v_max<=tmp(j)<=v_max];
    end
    
    %%%%%%%%%%%%%%% INTERVAL 1
    tmp=[viM1 a b]*Mvel_bs2basis2;
    tmp=tmp(:);
    for j=1:length(tmp)
        constraints=[constraints -v_max<=tmp(j)<=v_max];
    end   

    %%%%%%%%%%%%%%% INTERVAL 2
%     tmp=[a b c]*Mvel_bs2basis2;
%     tmp=tmp(:);
%     for j=1:length(tmp)
%         constraints=[constraints -v_max<=tmp(j)<=v_max];
%     end   
%  
%     %%%%%%%%%%%%%%% INTERVAL 3
%     V_MV3=sdpvar(3,3,'full');
%     constraints=[constraints V_MV3==[b c d]*Mvel_bs2basis3];
%     
%     tmp=V_MV3(:);
%     for j=1:length(tmp)
%         constraints=[constraints -v_max<=tmp(j)<=v_max];
%     end   
%  
%     %%%%%%%%%%%%% INTERVAL 4
%     V_MV4=sdpvar(3,3,'full');
%     constraints=[constraints V_MV4==[c d e]*Mvel_bs2basis4];
%     
%     tmp=V_MV4(:);
%     for j=1:length(tmp)
%         constraints=[constraints -v_max<=tmp(j)<=v_max];
%     end   
% 
%     %%%%%%%%%%%%% INTERVAL -2
%     V_MV5=sdpvar(3,3,'full');
%     constraints=[constraints V_MV5==[d e zeros(3,1)]*Mvel_bs2basis5];
%     
%     tmp=V_MV5(:);
%     for j=1:length(tmp)
%         constraints=[constraints -v_max<=tmp(j)<=v_max];
%     end   
%  
%     %%%%%%%%%%%%% INTERVAL -1
%     V_MV6=sdpvar(3,3,'full');
%     constraints=[constraints V_MV6==[e zeros(3,1) zeros(3,1)]*Mvel_bs2basis6];
%     
%     tmp=V_MV6(:);
%     for j=1:length(tmp)
%         constraints=[constraints -v_max<=tmp(j)<=v_max];
%     end  
    
    %%%%%%%%%%%%%%%
    clc;
    obj=-sum(a);
    result=optimize(constraints,obj,sdpsettings('usex0',0,'solver','fmincon','showprogress',0,'verbose',0,'debug',0 ));   
    value(a(1))
    
    obj=sum(a);
    result=optimize(constraints,obj,sdpsettings('usex0',0,'solver','fmincon','showprogress',0,'verbose',0,'debug',0 ));   
    value(a(1))
%     [ value(a)   value(b)   value(c)   value(d)  value(e)]
    
    
%     Vbs=[Vbs_first_block value(a)];
%     
%     Vbs*Mvel_bs2basis1;

%%
Mvel_basis2bs=inv(Mvel_bs2basis3);

vel_bs=[];
vel_mv=[];
increm=1.0;
for vx=-v_max:increm:v_max
    for vy=-v_max:increm:v_max
        for vz=-v_max:increm:v_max
            tmp=[vx vy vz]*Mvel_basis2bs;
            vel_bs=[vel_bs tmp'];
            vel_mv=[vel_mv [vx vy vz]'];
        end
    end
end

figure;
subplot(1,2,1);
scatter3(vel_bs(1,:),vel_bs(2,:),vel_bs(3,:),'y'); hold on;
scatter3(vel_mv(1,:),vel_mv(2,:),vel_mv(3,:),'g'); axis equal;
xlabel('viM2');
ylabel('viM1');
zlabel('vi');



Mvel_bs2basis=computeMatrixForClampedUniformBSpline(2,3,"01")*inv(getA_Be(2,"01"))
Mvel_basis2bs=inv(Mvel_bs2basis);

vel_bs=[];
vel_be=[];
increm=1.0;
for vx=-v_max:increm:v_max
    for vy=-v_max:increm:v_max
        for vz=-v_max:increm:v_max
            tmp=[vx vy vz]*Mvel_basis2bs;
            vel_bs=[vel_bs tmp'];
            vel_be=[vel_be [vx vy vz]'];
        end
    end
end

subplot(1,2,2);
scatter3(vel_bs(1,:),vel_bs(2,:),vel_bs(3,:),'y'); hold on;
scatter3(vel_be(1,:),vel_be(2,:),vel_be(3,:),'b'); axis equal;
xlabel('viM2');
ylabel('viM1');
zlabel('vi');
