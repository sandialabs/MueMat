% Drop entries from a matrix.
%
% This dropping function can be used as a PreDropFunc or
% PostDropFunc in CoalesceDropFactory.
%
% Foreach a_ij, we keep a_ij iff
%  (a_ij^2 / abs(a_ii*a_jj)) > tol^2 or
%  (a_ji^2 / abs(a_ii*a_jj)) > tol^2
%
% See also: SAdropExample
%
function [matrix] = SAdrop(matrix, MyLevel, tol)
  n = size(matrix,2);

  invDiag = abs(1./diag(matrix));
  invDiag = sparse(1:n,1:n,invDiag);

  pattern = matrix;
  pattern = pattern.^2;
  pattern = pattern*invDiag;     % col scaling
  pattern = (pattern'*invDiag)'; % row scaling

  toltol  = tol*tol;
  pattern = (pattern>toltol);

  pattern = spones(pattern+pattern');

  matrix = matrix .* pattern;

end