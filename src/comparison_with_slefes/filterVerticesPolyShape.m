function vertices=filterVerticesPolyShape(poly)

    vertices=poly.Vertices;
    vertices(any(isnan(vertices), 2), :) = []; %remove the rows with nnan
    vertices=uniquetol(vertices,'ByRows',true); %See https://www.mathworks.com/help/matlab/ref/uniquetol.html

end