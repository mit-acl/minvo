% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

% This file checks whether the centroid of each facet of the simplex is
% contained in the . 
% If this condition is true, the vector distance_to_each_centroid should have very small values
% after running this file
% Currently the condition is satisfied for degree <=3 (for higher, I think
% I may be running into numerical errors and tolerances from the solver)

close all; clc; clear;

addpath(genpath('./solutions'));

degree=3;
% V=rand(degree,degree+1);

V=[zeros(degree,1) eye(degree)]; %Standard simplex


interv=[-1,1];

P=V*getA_MV(degree,interv)

tan_points=getTangencyPoints(degree,interv);

tan_points_and_extrema=[tan_points; min(interv); max(interv)];

syms t real
T=[];
for i=0:(degree)
   T=[t^i ;T];
end


centroids=[]; %each column will contain a centroid of a facet
for i=1:size(V,2)
    V_tmp=V;
    V_tmp(:,i)=[]; %delete that vertex
    centroids=[centroids mean(V_tmp,2) ]; %And compute the mean of the rest of the vertexes
end

distance_to_each_centroid=realmax*ones(1,size(centroids,2));


[A rootsA]=getA_MV(degree, interv);
tan_points_and_extrema=cell2mat(rootsA)';

for i=1:size(tan_points_and_extrema,1)
    disp(i/size(tan_points_and_extrema,1)) %Percentage done
    for ii=1:size(tan_points_and_extrema,1)
        a=subs(P*T,t,tan_points_and_extrema(i)); %one extremum of the segment
        b=subs(P*T,t,tan_points_and_extrema(ii)); %another extremum of the segment

        for j=1:size(centroids,2) %And now compute the distance from that segment to each of the centroids
            
           centroid=centroids(:,j);
           distance= vpa(minimumDistancePointToSegment(a,b,centroid));
           distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance);
           
%            distance_other=realmax;
%            for tt=0:0.01:1
%                point_in_segment=a + tt * (b - a);
%                distance_other=min(distance_other,norm(point_in_segment-centroid));
%            end
%            vpa(distance_other,4)
%            
%            distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance_other);
           
        end

    end

end

distance_to_each_centroid


%%%%%%%%%%%%%%%%%%%%%

% for extrema=[-1,1]
%     a=subs(P*T,t,extrema); %one extremum of the segment
%     for tt=-1:0.001:1
%         tt
%         
%         b=subs(P*T,t,tt); %another extremum of the segment
%         
%         for j=1:size(centroids,2) %And now compute the distance from that segment to each of the centroids
%            distance= vpa(minimumDistancePointToSegment(a,b,centroids(:,j)));
%            distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance);
%         end
%         
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%



function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

%https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
function min_distance=minimumDistancePointToSegment(v,w,p) %   v|--------|w   and p is the point
    l2=norm(v-w)^2;
    t=max(0.0,min(1.0,(p-v)'*(w-v)/l2)); % We clamp t from [0,1] to handle points outside the segment vw.
    projection = v + t * (w - v);  % Projection falls on the segment
    min_distance=norm(p-projection);
end
