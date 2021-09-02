%Author: Jesus Tordesillas, jtorde@mit.edu, September 2021

%This function (and hence "computeSlefeScalar" as well) has been tested (to check for correctness) against the result
%produced by the "uniexample.c" function of http://www.cise.ufl.edu/research/SurfLab/download/SubLiME.tar.gz

function breakpoints=computeSlefe(P, num_seg, interv)

deg=size(P,2)-1; %degree
dim=size(P,1); %number of coordinates

num_of_breakpts=(num_seg+1);

p_down=[]; p_up=[];

for i=1:dim %For each of the coordinates
     [t_break_points, p_down_i, p_up_i]=computeSlefeScalar(P(i,:), deg, num_seg, interv);
     p_down=[p_down; p_down_i];
     p_up=[p_up; p_up_i];
end

%Now, p_down is a dim x num_of_breakpoints matrix  
%     p_up   is a dim x num_of_breakpoints matrix  

for j=1:num_of_breakpts
    
       %Now let's create a grid with all the combinations
       for i=1:dim %For each of the coordinates
             seg_j_coord{i}=  [p_down(i,j); %all the coordinates for the breakpoint j
                              p_up(i,j) ];
       end
       grid = cell(1,dim);
       [grid{:}] = ndgrid(seg_j_coord{:});

       %Let's save them as a dim x num_of_vertices matrix
       vertices=[];
       for i=1:dim
           vertices=[vertices; grid{i}(:)'];
       end
       %Now vertices is a dim x num_of_vertices matrix

       %Let's store it 
       breakpoints{j}.vertices=vertices;
       breakpoints{j}.t=t_break_points(j);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: In 2D, the code of this function would simply be

% [t_break_points, p_down_x, p_up_x]=computeSleveScalar(P(1,:), deg, num_seg, interv);
% [t_break_points, p_down_y, p_up_y]=computeSleveScalar(P(2,:), deg, num_seg, interv);
% 
% for i=1:(seg+1)
%    breakpoints{i}.vertices=[p_down_x(i)      p_up_x(i)    p_up_x(i)      p_down_x(i);       
%                            p_down_y(i)      p_up_y(i)    p_down_y(i)    p_up_y(i)];
%                        
%    breakpoints{i}.t=t_break_points(i);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
end