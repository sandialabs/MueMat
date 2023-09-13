function [A,dropcnt] = drop(X,tol,output,replaceVal)
%
% syntax:  [A,dropcnt] = drop(X,tol,output,replaceVal)
%
% Drop entries from X whose absolute values are less than or equal to tol.
% Returns X with appropriate entries dropped.
% Prints information if 'output' is nonzero.
% Note: this function does not reduce the memory requirement for A, i.e., the
%       new nonzeros are still stored in the matrix.  This makes this function
%       much faster then it would otherwise be.

[m,n] = size(X);
if (~varexist('output') || (varexist('output') && isempty(output)))
  output = 1;
end
if (~varexist('tol') || (varexist('tol') && isempty(tol)))
  tol = eps;
end

if (~varexist('replaceVal') || (varexist('replaceVal') && isempty(replaceVal)))
  replaceVal = 0;
end


if (nargout > 0)
  [aaa,bbb,ccc] = find( X );
  dropcnt = 0;
  for ii=1:length(aaa)
    if (abs(ccc(ii)) < tol) ccc(ii) = replaceVal; dropcnt = dropcnt+1; end
  end
  A = sparse(aaa,bbb,ccc,m,n);
else
  dropcnt = nnz(X);
  [aaa] = find( abs(X) > tol );
  dropcnt = dropcnt - length(aaa);
end

if (output)
  fprintf(1,'dropped %d of %d entries (%d%%), tol = %5.4e\n',dropcnt,nnz(X),round(100 * dropcnt/nnz(X)),tol);
end
