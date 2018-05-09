function [ result ] = ternaryOp( condition, resultIfTrue, resultIfFalse )
% Equivalent of ternary/conditional operator in C#
%   condition: your boolean condition
%   resultIfTrue:  the output if your condition is true
%   resultIfFalse: the output if your condition is false
result = resultIfFalse;
if condition
    result = resultIfTrue;
end
end

