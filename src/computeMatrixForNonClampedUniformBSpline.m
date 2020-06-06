function Abs=computeMatrixForNonClampedUniformBSpline(deg,interval)

    knots=0:15;
    segment_index=0; %The result is the same for all the segments
    computeMatrixForAnyBSpline(deg,segment_index,knots,interval) %Not

end
