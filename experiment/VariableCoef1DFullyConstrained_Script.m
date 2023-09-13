% See also: VariableCoef1DFullyConstrained.m

% Problem size
n = 81;

% Rho coefficients
rho = ones(n+1,1);
epsi = 1e-5;
i = floor(n/3);   rho(i) = epsi; rho(i+1) = epsi; % first jump
i = floor(2*n/3); rho(i) = epsi; rho(i+1) = epsi; % second jump

P = VariableCoef1DFullyConstrained(n, rho);

plot(P, 'LineWidth',3);