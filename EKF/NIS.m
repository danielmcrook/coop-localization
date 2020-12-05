function epyk = NIS(eyk,H,Pmk,R)
%NIS Summary of this function goes here
%   Detailed explanation goes here

sk = H*Pmk*H' + R;
Sk = 0.5*(sk + sk');


epyk = eyk'*inv(Sk)*eyk;

end

