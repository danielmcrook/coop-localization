function xk = correct(xk)
%CORRECT Summary of this function goes here
%   Detailed explanation goes here


% correct rad data
if xk(3) > pi
    xk(3) = xk(3) - 2*pi;
elseif xk(3) < -pi
    xk(3) = xk(3) + 2*pi;
end

if xk(6) > pi
    xk(6) = xk(6) - 2*pi;
elseif xk(6) < -pi
    xk(6) = xk(6) + 2*pi;
end




end

