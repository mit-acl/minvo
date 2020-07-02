close all; clc;
determs=[];
for n=1:7
    determs=[determs abs(det(getSolutionA(n,"m11")))]
end

plot(determs)

vpa(determs(3)/determs(2))

vpa(determs(4)/determs(2))

vpa(determs(5)/determs(2))

vpa(determs(6)/determs(2))

vpa(determs(7)/determs(2))