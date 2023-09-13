% Build the matrix of a 1D Laplace problem
%
% The matrix represents a finite difference discretization of the
% poisson problem on a uniform grid and with Dirichlet boundary
% conditions.
%
% Using the stencil notation, the matrix is of the form [-1 2 -1].
%
% Parameter 'n' is the number of grid point.
% The function return a scalar sparse matrix of size (n x n).
%
% Example: BuildLaplace1D(4) returns the following matrix:
%
%      2    -1     0     0
%     -1     2    -1     0
%      0    -1     2    -1
%      0     0    -1     2
%
% See also: BuildLaplace1DBlk, BuildLaplace2D

function [Amat] = BuildLaplace1D(n)

      fprintf('\n*****************************************************************************\n\n');
      fprintf('Building %d x %d scalar matrix\n',n,n);

      newMap = Map(n, 1);

      MatrixData = spalloc(n,n,3*n);
      for k=1:n
      MatrixData(k,k) = 2;
      if k ~= 1, MatrixData(k,k-1)= -1; end
      if k ~= n, MatrixData(k,k+1)= -1; end
      end

      %Amat.RowMap     = newMap;
      %Amat.ColMap     = newMap;
      %Amat.Apply      = @MatlabApply;
      %Amat.MatrixData = MatrixData;
      Amat = Operator(MatrixData,newMap,newMap,@MatlabApply, 'Laplace1D');
