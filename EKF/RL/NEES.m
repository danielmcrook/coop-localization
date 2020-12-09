function epxk = NEES(exk,Ppk)
%NEES Summary of this function goes here
%   Detailed explanation goes here

epxk = exk'*inv(Ppk)*exk;

end

