% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

% This file checks whether the centroid of each facet of the simplex belongs to the convex hull of the curve. 
% If this condition is true, the values of resnorm printed should be very small 
% What this file does is: 
    %1) we compute the convex hull of the curve (a polyhedron)
    %2) we compute the matrices A, b such that the polyhedron is Ax<=b
    %3) we compute the distance between that polyhedron and the centroids of the simplex by solving a linear least-squares problem

%More info available at https://www.mathworks.com/matlabcentral/answers/108655-convex-hull-higher-dimension
    
close all; clc; clear;  addpath(genpath('./../solutions')); addpath(genpath('./../utils'));

degree=4; interv=[-1,1];

V=[zeros(degree,1) eye(degree)]; %Standard simplex
P=V*getA_MV(degree,interv);

syms t real
T=getT(degree,t);

tt=min(interv):0.05:max(interv);
samples_curve=double(subs(P*T,t,tt));


centroids=[]; %each column will contain a centroid of a facet
for i=1:size(V,2)
    V_tmp=V;
    V_tmp(:,i)=[]; %delete that vertex
    centroids=[centroids mean(V_tmp,2) ]; %And compute the mean of the rest of the vertexes
end

%See https://www.mathworks.com/matlabcentral/answers/108655-convex-hull-higher-dimension

[A,b,Aeq,beq]=vert2lcon(samples_curve', 1e-10); %Ax<=b is the inside of the polyhedron

options = optimoptions('lsqlin','Algorithm','interior-point','Display','off');

for j=1:size(centroids,2)
%     disp("Starting");
    centroid=centroids(:,j);
    [x,resnorm,residual,exitflag,output,lambda]=lsqlin(eye(degree),centroid,A,b,Aeq,beq,[],[],[],options);
    fprintf('Squared distance to centroid %d = %8.8f\n',j,resnorm)
%     assert(resnorm<1e-4)
end

%%
%%%%%%%%%%%%%%%%%%%%%

% % OLD: 
% % This file checks whether the centroid of each facet of the simplex is
% % contained in the . 
% % If this condition is true, the vector distance_to_each_centroid should have very small values
% % after running this file
% % Currently the condition is satisfied for degree <=3 (for higher, I think checking the segments that join two extrema is not enough)
% 
% close all; clc; clear; addpath(genpath('./../solutions'));addpath(genpath('./../utils'));
% 
% degree=4;
% V=[zeros(degree,1) eye(degree)]; %Standard simplex
% interv=[-1,1];
% P=V*getA_MV(degree,interv);
% 
% syms t real
% T=getT(degree,t);
% PT=P*T;
% 
% 
% centroids=[]; %each column will contain a centroid of a facet
% for i=1:size(V,2)
%     V_tmp=V;
%     V_tmp(:,i)=[]; %delete that vertex
%     centroids=[centroids mean(V_tmp,2) ]; %And compute the mean of the rest of the vertexes
% end
% 
% distance_to_each_centroid=realmax*ones(1,size(centroids,2));
% 
% 
% [A rootsA]=getA_MV(degree, interv);
% tan_points_and_extrema=cell2mat(rootsA)';
% 
% for i=1:size(tan_points_and_extrema,1)
%     disp(i/size(tan_points_and_extrema,1)) %Percentage done
%     for ii=1:size(tan_points_and_extrema,1)
%         a=subs(PT,t,tan_points_and_extrema(i)); %one extremum of the segment
%         b=subs(PT,t,tan_points_and_extrema(ii)); %another extremum of the segment
% 
%         for j=1:size(centroids,2) %And now compute the distance from that segment to each of the centroids
%             
%            centroid=centroids(:,j);
%            distance= vpa(minimumDistancePointToSegment(a,b,centroid));
%            distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance);
%            
% %            distance_other=realmax;
% %            for tt=0:0.01:1
% %                point_in_segment=a + tt * (b - a);
% %                distance_other=min(distance_other,norm(point_in_segment-centroid));
% %            end
% %            vpa(distance_other,4)
% %            
% %            distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance_other);
%            
%         end
% 
%     end
% 
% end
% 
% distance_to_each_centroid


%%%%%%%%%%%%%%%%%%%%%
%%
% distance_to_each_centroid=realmax*ones(1,size(centroids,2));
% for extrema=[-1,1]
% %     a=subs(P*T,t,extrema); %one extremum of the segment
%     a=evalPT(P,extrema);
% 
% %     a-subs(P*T,t,extrema)
% %     assert(a-)
%     
%     for tt=-1:0.001:1
%         tt
%         b=evalPT(P,tt);
% %         b=subs(P*T,t,tt); %another extremum of the segment
%         
%         for j=1:size(centroids,2) %And now compute the distance from that segment to each of the centroids
%            distance= vpa(minimumDistancePointToSegment(a,b,centroids(:,j)));
%            distance_to_each_centroid(j)=min(distance_to_each_centroid(j), distance);
%         end
%         
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%

function result=evalPT(P,t)
    result=[];
    for i=1:size(P,1)
        result=[result; polyval(P(i,:),t)];
    end
end


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
