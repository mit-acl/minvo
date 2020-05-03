function polys=generateAdhocPolynomials(degree)

    polys=[];

    %Insert half of the polynomials

    deg_is_even = (rem(deg, 2) == 0);

    if(deg_is_even==0)

        B=sym('B',[((deg+1)/2),1],'real');
        R=sym('R',[(deg+1)/2,(deg-1)/2],'real');

        for i=1:((deg+1)/2)
            pol=-B(i)*(t-1);
            for j=1:(deg-1)/2
                pol=pol*((t-R(i,j))^2);
            end
            pol
            polys=[polys;pol];
        end
    else %Deg is even
        %WORK IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!

        B=sym('B',[(deg/2 +1),1],'real');
        R=sym('R',[deg/2 +1,(deg-1)/2],'real');

        for i=1:2:((deg)/2)
            pol=-B(i);
            for j=1:(deg)/2
                pol=pol*((t-R(i,j))^2);
            end
            polys=[polys;pol];
        end  

       for i=2:2:((deg)/2)
           pol=-B(i)*(t-1)*(t+1);
            for j=1:(deg-2)/2
                pol=pol*((t-R(i,j))^2);
            end
            polys=[polys;pol];
       end  

    end

%Insert the other half
polys=[polys;subs(polys,t,-t)]
    
%Solve now LP

end