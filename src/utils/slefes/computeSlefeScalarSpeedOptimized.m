%Author: Jesus Tordesillas, jtorde@mit.edu, December 2021
% This function is the same as computeSlefeScalar, but it is optimized for speed by avoiding the use of for/if statements and using Matlab's matrix operations

function [t_break_points, p_down, p_up, comp_time]= computeSlefeScalarSpeedOptimized(P, deg, num_seg, interv)

    %SECTION COMPUTED OFFLINE (i.e., NOT taken into account for comp time)
    %%%%%%%%%%%%%%%%%%%%%
    currentFile = mfilename( 'fullpath' );
    [pathstr,~,~] = fileparts( currentFile );
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
    all_ajdm=[];
    for j=1:(deg-1)    
        tmp=(j-1)*2*(num_break_points)+1;
        a_jdm_up=   data( tmp:(tmp+num_break_points-1) )';
        a_jdm_down= data( (tmp+num_break_points):(tmp+2*num_break_points-1) )';
        a_jdm=[data( tmp:(tmp+num_break_points-1) )     data( (tmp+num_break_points):(tmp+2*num_break_points-1) )];
        all_ajdm=[all_ajdm, a_jdm];
        
    end
    k=numel(t_break_points);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic;
    %%%%%%%%%%%%%%%%%% COMPUTATION TIME Starts
    all_second_diff=V -2 *V* diag(ones(numel(V)-1,1),-1) + V* diag(ones(numel(V)-2,1),-2);
    maxmin=[max(all_second_diff(1:deg-1),0.0);min(all_second_diff(1:deg-1),0.0)];
    l=((V(end)-V(1))/(t_break_points(end)-t_break_points(1)))*(t_break_points'-t_break_points(1)*ones(k,1)) + V(1)*ones(k,1); 
    p=[l l] + [all_ajdm*maxmin(:)   all_ajdm*reshape(flip(maxmin),[],1)];
    %%%%%%%%%%%%%%%%%% COMPUTATION TIME Ends
    comp_time=toc;

    p_down=p(:,2)';
    p_up=p(:,1)';

end