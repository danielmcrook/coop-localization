function epyk = NIS(eyk,H,Pmk,R)
%NIS Summary of this function goes here
%   Detailed explanation goes here

% cent = inv(chol(R,'lower'));
% cent = inv(chol(R,'lower'));

epyk = eyk'*inv(H*Pmk*H' + R)*eyk;

% for k=1:size(eyk,2)
%     
%     epyk(k) = eyk(:,k)'*cent*eyk(:,k);
%     
%     
%     
% end

end

