%Author: Jesus Tordesillas, jtorde@mit.edu, September 2021

%This code follows closely Sec. 3.3 of the paper "Threading Splines through 3D channels"

%Other related papers:
   % Efficient Pixel-Accurate Rendering of Curved Surfaces
   % SLEVEs for planar spline curves
   % Efficient Pixel-Accurate Rendering of Curved Surfaces
   % Package SubLiME (subdividable linear maximum-norm enclosure) package, downloadable from   http://www.cise.ufl.edu/research/SurfLab/download/SubLiME.tar.gz
   % https://vimeo.com/297544352#t=153s
   
%This function has been tested (to check for correctness) against the result pproduced by the "uniexample.c" function of http://www.cise.ufl.edu/research/SurfLab/download/SubLiME.tar.gz

   
function [t_break_points, p_down, p_up]= computeSlefeScalar(P, deg, num_seg, interv)


    %%% See  %See https://www.mathworks.com/matlabcentral/answers/250997-how-to-use-relative-path-to-use-matlab-file-in-another-computer
    currentFile = mfilename( 'fullpath' );
    [pathstr,~,~] = fileparts( currentFile );
    %%%
    
    filename=[pathstr,'/thirdparty/SubLiME/range/unirange-',num2str(deg),'_',num2str(num_seg),'.asc'];
    
    if ~isfile(filename)
        error("These values haven't been tabulated yet. Change the deg or # of segments")
    end
    
    data = readmatrix(filename,detectImportOptions(filename, 'FileType', 'text'));

    num_break_points=num_seg+1;
    t_break_points=linspace(min(interv),max(interv),num_break_points);

    assert(size(P,1)==1);

    %Compute the Bezier control points
    A_Be=getA_Be(deg,interv); 
    V=P*inv(A_Be);
    
%     l=(1-t_break_points)*V(1)+t_break_points*V(end); 
        
    %auxiliary variables
    x1=t_break_points(1);  x2=t_break_points(end);
    y1=V(1);               y2=V(end);
    
    l=((y2-y1)/(x2-x1))*(t_break_points-x1)  + y1;  %Linear interpolation between (x1,y1) and (x2,y2) (i.e., between the first and last Bezier control points)
                                                    %See paragraph above Eq. 1 in Sec. 3.3 of paper "Threading Splines through 3D channels"
    p_up=l;
    p_down=l;
    
    for j=1:(deg-1)    
        
        %auxiliary variables
        tmp=(j-1)*2*(num_break_points)+1;
        second_diff=(V(j)-2*V(j+1)+V(j+2)); %See Equation below Fig. 6 in Sec. 3.3 of paper "Threading Splines through 3D channels"
                                            %Note that Matlab uses 1-indexing
        
        
        a_jdm_up=   data( tmp:(tmp+num_break_points-1) )';
        a_jdm_down= data( (tmp+num_break_points):(tmp+2*num_break_points-1) )';
        
        %See equations below Eq. 1 in Sec. 3.3 of paper "Threading Splines through 3D channels"
        p_up=   p_up+ max(0,second_diff)*a_jdm_up + min(0,second_diff)*a_jdm_down; 
        p_down= p_down+ min(0,second_diff)*a_jdm_up + max(0,second_diff)*a_jdm_down; 
    end

end


%%
% for breakpoint_index=1:num_break_points
%     
%     %upper
%     t=t_break_points(breakpoint_index); %time at which that breakpoint happen
    
%     t=t_break_points;
    
% [t_break_points, p_down_x, p_up_x]=computeSleveScalar(V(1,:), deg, seg)

% 
% figure; hold on
% plot(t_break_points,V(1,:))
% plot(t_break_points,p_up_x)
% plot(t_break_points,p_down_x)
% A_Be=getA_Be(deg,[0,1]); syms t real
% P=V(1,:)*A_Be; T=getT(deg,t);
% fplot(P*T,[0,1]);
% end

% P=Q'
% [k,av] = convhull(P);
% plot(P(:,1),P(:,2),'*')
% hold on
% plot(P(k,1),P(k,2))
% patch(Q(1,:),Q(2,:),'red')

%     [k,av] = convhull(tmp');
%     plot(tmp(1,:),tmp(2,:),'*')
%     hold on
%     plot(tmp(1,k),tmp(2,k))
% plot(t_break_points,V)
% plot(t_break_points,p_up)
% plot(t_break_points,p_down)

% plot(p_down_x,p_down_y,'o')
% plot(p_up_x,p_up_y,'o')
% plot(p_up_x,p_down_y,'o')
% plot(p_down_x,p_up_y,'o')

% plot(t_break_points,p_down)