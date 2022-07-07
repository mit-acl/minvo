function [vol,surf] = objvol(obj)
% This function uses the "obj" structure created by the function readwObj.
% The calculation of volume and surface area of the 3d follows the
% algorithm presented in "Graphics Gems II" by James Arvo (ed.) page 170 -
% 171. The vector Qj was replaced by every face's centroid position vector.
% The formula for centroid position vector can be found in Wolfram
% Mathworld's "Geometric Centroid" article.
verts = obj.v; %obtain vertices' coordinates
faces = obj.f.v; %obtain face descriptor
area=zeros(size(faces,1),1);
centroid=zeros(size(faces,1),3);
normvec=zeros(size(faces,1),3);
for k=1:size(faces,1)
    sz=find(~isnan(faces(k,:)),1,'last');
    vectsid=cell(sz,sz-1);
    i=1;
    cntr=[0,0,0];
    while i<=sz-1
        if i==1;
            for j=1:sz
                vectsid{j,i}=verts(faces(k,j),:);
                cntr=cntr+verts(faces(k,j),:);
            end
        else
            for j=1:sz-i+1
                vectsid{j,i}=vectsid{j+1,i-1}-vectsid{1,i-1};
            end
        end
        i=i+1;
    end
    a=[0,0,0];
    for j=1:sz
        if j==sz
            a=a+cross(vectsid{1,1},vectsid{j,1});
        else
            a=a+cross(vectsid{j+1,1},vectsid{j,1});
        end
    end
    normvec(k,:)=cross(vectsid{1,sz-1},vectsid{2,sz-1});
    if all(normvec(k,:)==0);
        %some faces' vertex only form line (or zero-area polygons)
        normvec(k,:)=normvec(k,:);
    else
        %we calculate the area of faces that form valid polygons only.
        normvec(k,:)=normvec(k,:)./norm(normvec(k,:));
        area(k,1)=0.5*abs(dot(normvec(k,:),a));
        centroid(k,:)=cntr./3;
    end
end
volcom=dot(centroid,normvec,2).*area;
vol=sum(volcom)/3;
surf=sum(area);
end