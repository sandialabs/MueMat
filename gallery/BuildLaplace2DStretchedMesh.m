% Build the matrix of a 2D Laplace problem with stretched meshes
%
% The matrix represents a finite difference discretization of the
% poisson problem with the 5-point operator on a stretched mesh.
%
% Using the stencil notation, the matrix is of the form:
%
%      0    -1        0
%   -epsi  2+2*epsi -epsi
%      0    -1        0
%
% Parameters:
% 'n' is the number of grid point on each direction.
% 'epsi' is the stretch mesh ratio.
%
% The function return a scalar sparse matrix of size (n^2 x n^2).
%
% Example: BuildLaplace2D(2,0.3) returns the following matrix:
%
%     3.4   -1.0   -0.7     0.0
%    -1.0    3.4    0.0    -0.7
%    -0.7    0.0    3.4    -1.0
%     0.0   -0.7   -1.0     3.4
%
% See also: BuildLaplace2D

function [Amat] = BuildLaplace2DStretchedMesh(n, epsi)
% N is the total number of points.
N=n^2;

north = -1;
south = -1;
east  = -epsi;
west  = -epsi;
center = -(north+south+east+west);

newMap = Map(N, 1);

fprintf('\n*****************************************************************************\n\n');
fprintf('Building anisotropic %d x %d scalar matrix on stretched mesh...',N,N);

MatrixData = spalloc(N,N,4*N);

i=1; j=1;
for k=1:N
    % Center point
    MatrixData(k,k) = center;
    % North point
    if i ~= 1,  MatrixData(k,k-n) = north; end
    % South point
    if i ~= n,  MatrixData(k,k+n) = south; end
    % East point
    if j ~= 1,  MatrixData(k,k-1) = east; end
    % West point
    if j ~= n,  MatrixData(k,k+1) = west; end

    if j ~= n
        j=j+1;
    else
        i=i+1; j=1;
    end
end
%Amat.RowMap     = newMap;
%Amat.ColMap     = newMap;
%Amat.Apply      = @MatlabApply;
%Amat.MatrixData = MatrixData;

Amat = Operator(MatrixData,newMap,newMap,@MatlabApply, 'Laplace2DStretchedMesh');

fprintf(' done!\n');
fprintf('\n*****************************************************************************\n');
