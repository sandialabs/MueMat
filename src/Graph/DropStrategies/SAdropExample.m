% Example to write your own drop function for CoalesceDropFactory.
%
% The dropping criteria is the same in SAdrop and SadropExample but we used
% here a local function to decide if an entry a_ij must be dropped.
%
% Note: This code is dramatically slower than SAdrop and is only provide
% as an example to write your own drop function for the
% CoalesceDropFactory using an element-by-element viewpoint.
%
% See also: SAdrop

function [matrix] = SAdropExample(matrix, tol)
  [Arows,Acols,Avals] = find(matrix);
  NN = length(Arows);

  KeepIt = zeros(NN,1);

  for i=1:NN
    KeepIt(i) = KeepElement(Arows(i), Acols(i), matrix, tol);
  end

  Keepers = find(KeepIt);
  Arows   = Arows(Keepers);
  Acols   = Acols(Keepers);
  Avals   = Avals(Keepers);

  [n, m] = size(matrix);
  matrix = sparse(Arows,Acols,Avals, n, m);
end

% Foreach a_ij, we keep a_ij iff
%  (a_ij^2 / (a_ii*a_jj)) > tol^2 or
%  (a_ji^2 / (a_ii*a_jj)) > tol^2
function [keep] = KeepElement(i, j, matrix, tol)
  tol = tol*tol;
  t   = tol*abs(matrix(i,i))*abs(matrix(j,j)); % col/row diagonal scaling

  keep = 0;
  if (matrix(i,j)^2 >=  t) || (matrix(j,i)^2 >=  t), keep = 1; end
end

% Some ways to optimize such codes in matlab:
% - consider rewriting using a matrix point-of-view (see SAdrop)
% else
% - use Aval(i) instead of matrix(ii,jj) in KeepElement. It is significantly faster
% - store diagonal entries in a temporary vector to perform the scaling:
%   add matrixDiag=abs(diag(matrix)) to the arguments of KeepElement
% - compute tol*tol outside of KeepElement
% - take advantage of the symmetry of the dropping criteria