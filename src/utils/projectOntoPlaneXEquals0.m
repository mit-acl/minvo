function result=projectOntoPlaneXEquals0(p,m)
    proj=simplify(p(2:end)/(1-p(1)));  %Here we are projecting onto the plane x=0, using as viewpoint the vertex [1 0 0 ..... 0]
    if(numel(p)>(m+1))
        numel(proj)
       result=projectOntoPlaneXEquals0(proj,m);
    else
       result=proj;
    end
end