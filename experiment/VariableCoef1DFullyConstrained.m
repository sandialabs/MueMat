% The goal of this experiment is to demonstrates that constraints may not capture sufficiently physics.
%
% This example construct a 1D example with jumps in material
% properties and compute the prolongator:
% [     1    ]  (ie: RAP = Schur Complement)
% [ -Aff\Afc ]
%
% This prolongator correspond to fully constrain the problem with:
% nullspace = [ (1) (x) ] (but it's easier to just compute -Aff\Afc)
%
% Problem:
% ( rho^(x,y) u_x )_x + ( rho^(x,y) u_y )_y
%
% Discretisation:
%
% rho_3(u_3-u_2) - rho_2(u_2-u_1)
%       -------          -------
%         h                 h        rho_3 u_3 + rho_2 u_1 - (rho_3 + rho_2) u_2
% ------------------------------- = ---------------------------------------------
%                h                                        h^2
%
%
%        u1      u2      u3      u4      u5      u6      u7      u8      u9
% - - - -|- - - -|- - - -|- - - -|- - - -|- - - -|- - - -|- - - -|- - - -|- - - -
%  rho1    rho2    rho3    rho4    rho5    rho6    rho7    rho8    rho9    rho10
%
%
% Notes:
% --------
% Hard coded aggregation creates aggregates of size 9.
% Roots points are at the middle of aggregates.
%
% Usage example: See VariableCoef1DFullyConstrained_Script.m
%
% See also: VariableCoef1DFullyConstrained_Script.m
%
%
function P = VariableCoef1DFullyConstrained(n,rho)

%% Build A
A = sparse(n,n);
for i=1:n-1
   A(i,i)   = rho(i) + rho(i+1);
   A(i,i+1) = -rho(i+1);
   A(i+1,i) = -rho(i+1);
end
A(n,n) = rho(n) + rho(n+1);

%% Aggregation
%% - big aggregates (of size 9)
%% - roots points at the middle of aggregates
cpoints = zeros(n,1); cpoints(5:9:end) = 1; cpoints = find(cpoints);
fpoints = ones(n,1);  fpoints(cpoints) = 0; fpoints = find(fpoints);
nAggregates = size(cpoints,1);

Aff = A(fpoints,fpoints);
% Acf = A(cpoints,fpoints); % unused
% Acc = A(cpoints,cpoints);  % unused
Afc = A(fpoints,cpoints);

%% Compute P
Pf = -Aff\Afc;
P = spalloc(n,nAggregates,1);
P(cpoints,1:nAggregates) = eye(nAggregates);
P(fpoints,1:nAggregates) = Pf;

end