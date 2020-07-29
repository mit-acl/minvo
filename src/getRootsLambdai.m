function roots_lambdai=getRootsLambdai(degree, interval)
%TODO set correct order in lambda_i
if(degree>7 || degree<=0)
    error('not implemented yet!')
end

if(degree==1)
    roots_lambdai{1}=[1.0];
    roots_lambdai{2}=[-1.0];
else
    sol=load(strcat('sols_formula/rootsLambdaiDeg',num2str(degree),'.mat'));
    roots_lambdai=sol.roots_lambda_solution;
end
    
if(interval~="m11")
    error('not implemented yet!')
end

end