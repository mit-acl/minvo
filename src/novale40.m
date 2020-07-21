     

deg=7;

integrals_vector=[];
for deg=3:2:11
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

       expand(prod)
%      integrals=prod;
%      for i=((deg-1)/2-1):-1:0
%           integrals= int( integrals  ,tt(i+2),-1,1);
%      end
%      integrals= int( integrals  ,tt(1),-1,1)
%      
%      
%      
% integrals_vector=[integrals_vector integrals];
     
end