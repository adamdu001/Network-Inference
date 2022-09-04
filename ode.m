%ODE funtion
function theta_dot = ode(t,theta,omega,K,N,A)
         theta_dot = omega + K/N*sum(A'.*sin(theta-theta'))';
end
% theta-theta' gives a N by N matrix
% for each col of theta-theta', it's essentially j-i for all j
% apply a sum on this matrix is to sum each col, gives 1 by N vector
% that's why we transpose it back to N by 1

