function obj = readwObj(fname)
% To use, put the *.obj files in the subfolders or same folder with this
% file, and navigate the MATLAB's current folder into the one containing
% them. Use the syntax obj=readwObj('fname'), the fname is the full path to
% your obj file. I designed this obj file reader to deal with TEXTURELESS
% Wavefront Object file created by DAVID Laserscanner's Shapefusion. This
% file was created in MATLAB 2016 with possible compatibility to MATLAB
% 2014.
% Author: Alutsyah Luthfian
% Date: March 2, 2018
% If the faces in the *.obj file were not uniformly made from similar
% amount of vertices (face A made from 5 vertices, face B made from 3
% vertices, and so on), the face matrix in the "obj" structure will be of
% n*m size, n is the amount of faces of the 3d shape and m is the largest
% amount of vertices needed to make those faces. Other faces which made
% from less vertices will have zeros in the columns succeeding the last
% vertex. To calculate the volume and surface area of the 3d shape, please
% use the objvol function included in this package.
fid=fopen(fname);
numv=0;
numn=0;
numf=0;
lenf=zeros(2,1);
while ~feof(fid)
    tline=fgetl(fid);
    ln = sscanf(tline,'%s',1);
    switch ln
        case 'v'
            numv=numv+1;
        case 'vn'
            numn=numn+1;
        case 'f'
            numf=numf+1;
            face=sscanf(tline(2:end),'%i//%i',[2 Inf]);
            if length(face)==1&&~isempty(strfind(tline(2:end),'/'))
                face=sscanf(tline(2:end),'%i/%*i/%i',[2 Inf]);
            else
                face=sscanf(tline(2:end),'%i',[1 Inf]);
            end
            lenf(2,1)=length(face);
            if lenf(2,1)>lenf(1,1)
                lenf(1,1)=lenf(2,1);
            end
    end
end
obj.v=zeros(numv,3);
obj.vn=zeros(numn,3);
obj.f.v=nan(numf,max(lenf));
obj.f.vn=nan(numf,max(lenf));
frewind(fid);
numv=1; numn=1; numf=1;
while ~feof(fid)
    tline=fgetl(fid);
    ln = sscanf(tline,'%s',1);
    switch ln
        case 'v'
            obj.v(numv,:)=sscanf(tline(2:end),'%f',[3 Inf])';
            numv=numv+1;
        case 'vn'
            obj.vn(numn,:)=sscanf(tline(3:end),'%f',[3 Inf])';
            numn=numn+1;
        case 'f'
            face=sscanf(tline(2:end),'%i//%i',[2 Inf]);
            if length(face)==1&&~isempty(strfind(tline(2:end),'/'))
                face=sscanf(tline(2:end),'%i/%*i/%i',[2 Inf]);
                obj.f.v(numf,1:length(face))=face(1,:);
                obj.f.vn(numf,1:length(face))=face(2,:);
                numf=numf+1;
            elseif length(face)==1
                face=sscanf(tline(2:end),'%i',[1 Inf]);
                obj.f.v(numf,1:length(face))=face;
                numf=numf+1;
            else
                obj.f.v(numf,1:length(face))=face(1,:);
                obj.f.vn(numf,1:length(face))=face(2,:);
                numf=numf+1;
            end
    end
end
fclose(fid);