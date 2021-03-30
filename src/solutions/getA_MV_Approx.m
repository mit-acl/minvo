function [A rootsA]=getA_MV_Approx(n,interv)

    rootsA=getRoots_MV_Approx(n,interv);

    %Now I need to compute the coefficients of the polynomial that satisfy A'1=ee

    ee=[zeros(1,n) 1]';

    %Option 1: Use symbolic toolbox
    % b=sym('b', [n+1 1],'real');
    % A=[];
    % for i=0:n
    %     A=[A; b(i+1)*poly(approx_mv_roots(indexes==i))];
    % end
    % b_solved=solve(sum(A)==[zeros(1,n) 1]);
    % A=double(subs(A,b_solved))


    %Option 2 (same as before, but solve "manually" the linear system avoiding
    % the symbolic toolbox)
    C=[];
    for i=1:(n+1)
        C=[C; poly(rootsA{i})];
    end

    b=linsolve(C',ee);  %Other options: b=C'\ee;  b=pinv(C')*ee;   TODO: this step can run into numerical issues if n is very large, because C is close to singular. 
    A=diag(b)*C;

    %Option 3
%     C=[];
%     for i=1:(n+1)
%         C=[C; poly(rootsA{i})];
%     end
%     b=sym('b', [n+1 1],'real');
%     tt=linspace(-1,1,n+1)'.^(n:-1:0);
%     b_solved=solve(sum(diag(b)*C*tt)==ones(1,n+1))
%     A=double(diag(subs(b,b_solved))*C);
    
    norm(sum(A)-ee')
%     assert(norm(sum(A)-ee')<1e-3);

end