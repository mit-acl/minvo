function result=convSym(a,b) %MATLAB conv() function doesn't work with symbolic variables

    syms t real
    nn=numel(a,2)-1;
    result=[];
    poly_multiplication=(a(:)'*getT(numel(a)-1,t))*(b(:)'*getT(numel(b)-1,t)); %Note that discrete convolution is the same as polynomial multiplication
    coefficients=coeffs(poly_multiplication,t,'All');
    coefficients=[zeros(1,2*nn-numel(coefficients)+1) coefficients];
    result=[result ;coefficients];
    
end