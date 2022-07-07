%See https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html
classdef MyAverageKeeper< handle

    properties
        average
        size_so_far
    end
    
    methods

        function obj = MyAverageKeeper() 
            obj.average=0.0;
            obj.size_so_far=0;

        end

        %https://stackoverflow.com/a/22999488/6057617
        function addToAverage(obj, value)
            obj.average=(obj.size_so_far * obj.average + value)/(obj.size_so_far+1);
            obj.size_so_far=obj.size_so_far+1;
        end
 
    end
end