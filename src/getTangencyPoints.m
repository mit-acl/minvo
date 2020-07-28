function tangency_points=getTangencyPoints(degree, interval)

if(degree>7 || degree<=0)
    error('not implemented yet!')
end

if(degree==1)
    tangency_points= []; 
else

sol=load(strcat('sols_formula/solutionTangencyPointsDeg',num2str(degree),'.mat'));
tangency_points=sol.tangencyPoints;

end
    
if(interval~="m11")
    error('not implemented yet!')
end

end