close all; clc; clear;
syms t

% global R B_solved determ

deg=4;
n=deg;

W=[];



%Insert half of the polynomials

deg_is_even = (rem(deg, 2) == 0);

if(deg_is_even==0)
    
    B=sym('B',[((deg+1)/2),1],'real');
    R=sym('R',[(deg+1)/2,(deg-1)/2],'real');
    index_B=1;
    
%     for i=1:((deg+1)/2)
%         pol=-B(i)*(t-1);
%         for j=1:(deg-1)/2
%             pol=pol*((t-R(i,j))^2);
%         end
%         pol
%         W=[W;pol];
%     end

    for i=1:2:deg
        pol=-B(index_B)*(t-1);
        for j=1:(deg-1)/2
            pol=pol*((t-R(index_B,j))^2);
        end
        pol
        W=[W;pol];
        index_B=index_B+1;
    end
    
    W=[W;subs(W,t,-t)] %Insert the other half
    
else %Deg is even

    n_2_is_even = (rem(n/2, 2) == 0);
    if(n_2_is_even)
        Ia=2:2:(n/2);
        Ib=1:2:(n/2); %n/2+1
    else
        Ia=2:2:(n/2); %n/2+1
        Ib=1:2:(n/2); 
    end
    
    B=[]; %sym('B',[(deg/2+1),1],'real');
    R=[];
%     R=sym('R',[deg/2+1,100],'real');
    
    for index=1:size(Ia,2)
        i=Ia(index);
        Bi=sym(strcat('B',num2str(i))); B=[B; Bi];
        pol=-Bi*(t+1)*(t-1);
        for j=1:(n-2)/2
            Rij=sym(strcat('R',num2str(i),num2str(j))); R=[R; Rij];
            pol=pol*((t-Rij)^2);
        end
        W=[W;pol];
        W=[W;subs(pol,t,-t)];
    end  
 
   for index=1:size(Ib,2)
       i=Ib(index);
       Bi=sym(strcat('B',num2str(i))); B=[B; Bi];
       pol=Bi;
       for j=1:(n/2)
            Rij=sym(strcat('R',num2str(i),num2str(j))); R=[R; Rij];
            pol=pol*((t-Rij)^2);
       end
       W=[W;pol];
       W=[W;subs(pol,t,-t)];

   end  
   
   %Now, for the labmda that is in the middle, force it to have symmetrical
   %roots
   i=n/2+1;
   if(n_2_is_even)
       Bi=sym(strcat('B',num2str(n/2+1))); B=[B; Bi];
       pol=B(i);
       for j=1:(n/4)
            Rij=sym(strcat('R',num2str(i),num2str(j))); R=[R; Rij];
            pol=pol*((t-Rij)^2)*((t+Rij)^2);
       end
       W=[W;pol];
   else
       Bi=sym(strcat('B',num2str(n/2+1))); B=[B; Bi];
       pol=B(i)*(t+1)*(t-1);
       for j=1:(n/4 - 2)
            Rij=sym(strcat('R',num2str(i),num2str(j)));R=[R; Rij];
            pol=pol*((t-Rij)^2)*((t+Rij)^2);
       end
       W=[W;pol];
   end
   
%    disp('W before')
%    W
%    disp('W after')
%    W=[W;subs(W(1:n/2,:),t,-t)] %Insert the other half
 
end


%Solve for the coefficients:
coeffic=coeffs(sum(W),t,'All'); %When using All, no flip!! (gives directly [a  b c d]

coeffic=[zeros(1,deg+1-length(coeffic)) coeffic];
%%
disp("Solving linear system")
solution=solve(coeffic==[zeros(1,deg) 1],B); %This is a linear system

% solve(coeffic(3:5)==[0 0 1],B);

disp("struct2array")
if(deg~=1)
    B_solved=(struct2array(solution)');
    B_solved=B_solved(:);
else
    B_solved=solution;
end

disp("Substituting")
W=subs(W,B,B_solved);

disp("Creating the A matrix")

%Create the A matrix:
A=[];
for i=1:length(W)
    tmp=coeffs(W(i),t,'All');
    A=[A; tmp];
end
disp("Computing the determinant")
%Compute the determinant
if(deg_is_even==0)
    determ=computeDetSmartly(A);
else
    determ=det(A)
end



% fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% 'lb',lb,'ub',ub,
lb=-ones(size(R));
ub=ones(size(R));

x0 = rand(size(R));
%%
opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',10000);
problem = createOptimProblem('fmincon','objective', ...
    @(R_x) getObj(R_x,R, determ),'x0',x0,'lb',-ones(size(R)),'ub',ones(size(R)),...
    'nonlcon',@(R_x) constraints(R_x,R, B_solved),'options',opts);



% Construct a GlobalSearch object
gs = GlobalSearch('Display','iter');
% Construct a MultiStart object based on our GlobalSearch attributes
ms = MultiStart('Display','iter','UseParallel',true);

disp('Running, it usually takes some time until the parpool starts');
[xgs,~,~,~,solsgs] = run(ms,problem,30); %8000

% [xms,~,~,~,solsgs] = run(ms,problem,100);

%%%%%%%%%%%%%%%%%%%%%%%%%
% ms = MultiStart(ms,'UseParallel',true);
% 
% try
%     demoOpenedPool = false;
%     % Create a parallel pool if one does not already exist
%     % (requires Parallel Computing Toolbox)
%     if max(size(gcp)) == 0 % if no pool
%         parpool
%         demoOpenedPool = true;
%     end
% catch ME
%     warning(message('globaloptim:globaloptimdemos:opticalInterferenceDemo:noPCT'));
% end
% 
% % Run the solver
% tic;
% [xms,~,~,~,solsgs] = run(ms,problem,100);
% toc
% xms

%%%%%%%%%%%%%%%%%%%%
%% 
B_solution=vpa(subs(B_solved,R,xgs));
A_solution=double(vpa(subs(A,R,xgs)));
det(A_solution)

A=A_solution;

syms t real
T=[];
for i=0:(deg)
   T=[t^i ;T];
end
interv=[-1,1];
figure;
fplot(A_solution*T,interv)

rootsA=[];
for i=1:size(A_solution,1)
    rootsA=[rootsA ; roots(A_solution(i,:))'];
end
rootsA=double(real(rootsA));
save(['solutionDeg' num2str(deg) '.mat'],'A','rootsA');

function result=getObj(R_x, R, determ)
    
%     result=det(subs(A,R,R_x));
    %global R determ
     result=subs(determ,R,R_x);
    result=-abs(double(result));
end

function [c,ceq] =constraints(R_x, R, B_solved)
    %global R B_solved
    c=-subs(B_solved(:),R,R_x); %  -B(i) <=0
%     c=[c ; subs(R(:),R,R_x) - 2*ones(size(R(:)))]; %R(i) <= 2
%     c=[c ; -subs(R(:),R,R_x) - 2*ones(size(R(:)))]; % R(i) >=-2 //// -2 - R(i) <=0
    c=double(c);
    ceq=[];
end

function result=computeDetSmartly(A)

%odd columns of A
oddcol=A(:,1:2:end);

%even columns 
evencol=A(:,2:2:end);

A_reorganized=[oddcol,evencol];

nrow=size(A,1);
ncol=size(A,2);

AA=A_reorganized(1:nrow/2,1:ncol/2);

BB=A_reorganized(1:nrow/2,ncol/2+1:end);

result=det(2*AA)*det(BB);

end

% function [c,ceq] = ellipseparabola(x)
% c(1) = (x(1)^2)/9 + (x(2)^2)/4 - 1;
% c(2) = x(1)^2 - x(2) - 1;
% ceq = [];
% end
