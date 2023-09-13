% Builds a distance Laplacian mimicing the stencil of the given
% matrix. NOTE: If you have multiple equations per node, this will
% not work.  Hand in the agglomerated matrix, if you please.
%
% [Amat] = BuildDistanceLaplacian(Bmat,coords,[dfunc])
%
% Input:
%  Bmat   - Matrix whose stencil to mimic
%  coords - Coordinates of nodes
%  [dfunc]- distance function (dist=dfunc(x,y)).  [Default=2-norm]
% Output:
%  Amat   - The distance Laplacian

function [Amat] = BuildDistanceLaplacian(Bmat,coords,varargin)
% Set distance measure
if(length(varargin)==0), dfunc=@(x_,y_)norm(x_-y_,2);
else dfunc=varargin{1};end

% Grab graph and calculate distance Laplacian (off-diagonals)
[II,JJ,VV]=find(Bmat.GetMatrixData());
[M,N]=size(Bmat.GetMatrixData());
VV=VV*0;
for K=1:length(VV),
  DIST=dfunc(coords(II(K),:),coords(JJ(K),:))^2;
  if(II(K)~=JJ(K)), VV(K)=1/DIST;end
end
% Get Diagonal Correct
A=sparse(II,JJ,VV,M,N);
A=A+spdiags(A*ones(N,1),0,M,N);

Amat = Operator(A,Bmat.GetRowMap(),Bmat.GetColMap(),@MatlabApply,'Distance Laplacian');
