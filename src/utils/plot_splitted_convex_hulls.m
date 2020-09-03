% /* ----------------------------------------------------------------------------
%  * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
%  * Massachusetts Institute of Technology
%  * All Rights Reserved
%  * Authors: Jesus Tordesillas, et al.
%  * See LICENSE file for the license information
%  * -------------------------------------------------------------------------- */

function volume=plot_splitted_convex_hulls(P,A,num_of_intervals,color,radius_sphere);

samples=[];
samples_t=linspace(-1,1,num_of_intervals+1);
all_vertexes=[];%Its columns are the vertexes
for i=1:(length(samples_t)-1)
    a=samples_t(i);
    b=samples_t(i+1);
    P_converted=convertPFromABtoCD(P,[a,b],[-1,1]);
    V=P_converted*inv(A);
    all_vertexes=[all_vertexes V];
    %plot_convex_hull(P_converted(1,:)',P_converted(2,:)',P_converted(3,:)',A,'b',0.0017);    
end

color_vertex=[.98 .45 .02];

for i=1:size(all_vertexes,2)
    s1=plotSphere(all_vertexes(:,i),radius_sphere, color_vertex);
end

axis equal
tmp=gca;
if (size(findobj(tmp.Children,'Type','Light'))<1) %If still no light in the subplot
 camlight %create light
end
lighting phong

 
x=all_vertexes(1,:); y=all_vertexes(2,:); z=all_vertexes(3,:);
[k1,volume] = convhull(x,y,z);
s2=trisurf(k1,x,y,z,'LineWidth',1,'FaceColor',color);
alpha(s2,0.1)

end