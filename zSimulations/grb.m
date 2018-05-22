function [ out ] = grb( e,timeVal,h )
% get the appropriate value from a predetermined distribution
% e = distribution whose length is the maximum number of values to check
% t = current point in time
% h = step size (per runge kutta - which actually means you'll use some h/2
% terms as well...so double the "true" step size)
%
% out = output value from the distribution

index = round((timeVal*2/h)+1);
out = e(index);

end

