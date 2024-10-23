function y = piecewiseLine2(x,ton,ts1,toff,ka,kb)
% PIECEWISELINE   A line made of five pieces
% that is not continuous.

y = zeros(size(x));
ts2=-(ka*(ton-ts1)-kb*toff)/kb;
% This example includes a for-loop and if statement
% purely for example purposes.
for i = 1:length(x)
    if x(i) < ts1
        y(i) = 0;
    elseif x(i)>ts1&&x(i)<ton
        y(i) = ka*(x(i)-ts1);
    elseif x(i)>=ton&&x(i)<=toff
        y(i) = ka*(ton-ts1);
    elseif x(i)>toff&&x(i)<ts2
        y(i) = ka*(ton-ts1)+kb*(x(i)-toff);
    elseif x(i)>=ts2
        y(i) = 0;
    end
end
end