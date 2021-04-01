function rootsA=getRoots_MV_Approx(n,interv)

   %Compute all the inner roots in interv
   %-------------------------------------
   
   %These values have been found in fit_MV_roots.m
   b=[0.27350297910476373264287985875853, 3.0385093464854278089148920116713, 0.47785675906936553314352522647823];%, 772.6680998615668158890912309289
   approx_roots=[];
   for ic=0:(n-1) %index of the cluster
       numElemClusterm1=(floor((n+(isOdd(ic) & isOdd(n)))/2)) -1; %Number of elements of the cluster -1
       for iic=0:numElemClusterm1 %index inside the cluster
            root_i=sin((b(1)*(iic-(numElemClusterm1)/2.0)    +  b(2)*(ic-(n-1)/2.0))./(n+b(3))); %This uses the interval [-1,1] 
            approx_roots=[approx_roots   root_i];
       end
   end
   
    %Compute the index of the polynomial of the basis (0,...,n) that has those roots. I.e.: Polynomial whose index is indexes(i) will have root approx_roots(i)
    %---------------------------------------------------------------------------------------------------------------------------------------------------------

    cluster{1}=num2cell(2:2:(n));
    for j=2:n %for each cluster
        cluster{j}={};
        for jj=1:numel(cluster{j-1})
            a=getNext(cluster{j-1}{jj},n); %Note that getNext was designed to match the results obtained for n=1...,7
            cluster{j}{end+1}=a;
        end
    end
    
    indexes=[];
    for j=1:n
        indexes=[indexes cell2mat(cluster{j})];
    end
    
    %%%%% Double roots
    approx_roots=[approx_roots approx_roots]; %duplicate because they are double roots    
    indexes=[indexes indexes]; %duplicate because they are double roots
    
    %%%%Roots at -1
    approx_roots=[approx_roots -ones(1,ceil(n/2))]; %roots at -1;
    indexes=[indexes  1:2:n]; %polys that have roots at -1;
    
    %%%%Roots at +1
    approx_roots=[approx_roots  ones(1,ceil(n/2))]; %roots at +1;   
    if(isEven(n))
        indexes=[indexes  1:2:n]; %roots at +1;
    else
        indexes=[indexes  0:2:n]; %roots at +1;
    end
    
    
    %Change the interval
    for j=1:length(approx_roots)
         approx_roots(j)=convertNumberFromABtoCD(approx_roots(j),[-1,1],interv);
    end
    
    for j=0:n
        rootsA{j+1}=approx_roots(indexes==j);
        rootsA{j+1}=sort(rootsA{j+1}); 
    end
    
end

%%If n is even, this generates this sequence: 0 -> 1 -> 2 -> 3 -> ... -> n -> 0 -> 1 -> 2 ->....  (all the elements or the list are scalars
%%If n is odd, this generates this sequence:  1 -> 2 -> 3 ->...-> n-1 ->[n,0] -> 1 -> 2  ->  .... (all are scalars except [n,0])
%Note also that getNext was designed to match the results obtained for n=1...,7
function result=getNext(i,n)
   if(isEven(n))
       result=rem(i+1,n+1);
   else
       if(i(1)==n)
           result=1;
       elseif(i==(n-1))
            result=[n,0];
       else
           result=rem(i+1,n+1);
       end
   end
end