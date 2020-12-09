function yk = NLmeas(x)
%NLMEAS Summary of this function goes here
%   Detailed explanation goes here

xig = x(1);
etag = x(2);
thetag = x(3);
xia = x(4);
etaa = x(5);
thetaa = x(6);

yk = [atan2((etaa-etag),(xia-xig))-thetag...
    sqrt((xig-xia)^2 + (etag-etaa)^2)...
    atan2((etag-etaa),(xig-xia))-thetaa...
    xia etaa]';

end

