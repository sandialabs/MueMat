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
function [matrix] = ModDrop(matrix, MyLevel, tol)
  n = size(matrix,2);
  
  invDiag = abs(1./diag(matrix));
  invDiag = sparse(1:n,1:n,invDiag);
  
  pattern = matrix;
NnzPerRow = ((spones(pattern)*ones(n,1)).^2.0);
AvgNzPerRow = sum(NnzPerRow)/n;
  pattern = pattern.^2;
  pattern = pattern*invDiag;     % col scaling
  pattern = (pattern'*invDiag)'; % row scaling
% find largest off-diag nonzero
%temp  = pattern - spdiags(diag(pattern),0,n,n);
%RowMax = max(abs(temp'))';
%tattern = spdiags(1./RowMax,0,n,n)*temp + speye(n,n);
  
  toltol  = tol*tol;
[aaa,bbb,ccc] = find(pattern);
ccc = (ccc .* NnzPerRow(aaa))/AvgNzPerRow;
inds = find( abs(ccc) > toltol);
pattern = sparse(aaa(inds),bbb(inds),ones(length(inds),1),n,n) + speye(n,n);

%  pattern = (tattern>toltol);
  
  pattern = spones(pattern+pattern');
  
  matrix = matrix .* pattern;

end
