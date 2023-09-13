%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function MAT=ssubsasgn(MAT,I,J,MAT2)
% Assigns MAT(I,J) = MAT2, for vectorized indices I and J.  Empty
% arguements are acceptable as replacements for the : which can't
% be used here.
%
% Why write this function?  Because somewhere in its dark innards,
% MATLAB's subsasgn function for sparse matrices makes everything
% dense, causing a *massive* use of memory.  This is, of course,
% documented absolutely nowhere.  Thanks, Mathworks!
%
% by Chris Siefert <csiefer@sandia.gov>
% Last Updated: 05/01/06
function MAT=ssubsasgn(MAT,I,J,MAT2)
[M,N]=size(MAT);
[K,L,V]=find(MAT2);
if(isempty(I)), MAT=MAT+sparse(K,J(L),V,M,N);
elseif(isempty(J)), MAT=MAT+sparse(I(K),L,V,M,N);
else MAT=MAT+sparse(I(K),J(L),V,M,N);
end
