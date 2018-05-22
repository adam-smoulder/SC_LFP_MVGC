function [ Q ] = qf( g, Qm0 )
h0 = -log(1+log(1+1/Qm0));
if g > -h0
    Q = Qm0*(1-exp(-(exp(g)-1)/Qm0));
else
    Q = -1;
end
end

