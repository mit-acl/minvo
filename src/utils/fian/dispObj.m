function [verts, faces]=dispObj(obj)
%
% function display_obj(obj)
%
% displays a TEXTURELESS obj structure
%
% INPUTS:  obj:     object data
%                   - obj.v:    vertices
%                   - obj.vt:   texture coordinates
%                   - obj.f.v:  face definition vertices
%                   - obj.f.vt: face definition texture
%
% Modified by Alutsyah Luthfian (2018)
% Original Author: Bernard Abayowa
% University of Dayton
% 6/16/08
cntr=mean(obj.v,1);
tval=zeros(size(obj.v,1),1);
for i=1:size(obj.v,1)
    tval(i,1)=0.0;%*norm(obj.v(i,:)-cntr);
end
% display object
% figure
p=patch('vertices',obj.v,'faces',obj.f.v,'FaceVertexCData', tval);
shading interp
colormap jet;
lighting phong;
camlight('right');
camproj('perspective');
view(45,45)
axis square;
axis off;
axis equal
axis tight;
%cameramenu
verts=get(p,'Vertices');
faces=get(p,'Faces');
