function [Amat] = BuildLaplace2D(n,periodic)
% Build the matrix of a 2D Laplace problem
%
% syntax [Amat] = BuildLaplace2D(n,periodic)
%
% This function is a rewriting of the matlab built-in
% function gallery('poisson',n).
%
% The matrix represents a finite difference discretization of the
% poisson problem with the 5-point operator on a uniform grid and
% with Dirichlet boundary conditions.
%
% Using the stencil notation, the matrix is of the form:
%
%      0    -1     0
%     -1     4    -1
%      0    -1     0
%
% Parameter 'n' is the number of grid point on each direction.
% The function return a scalar sparse matrix of size (n^2 x n^2).
%
% Example: BuildLaplace2D(2) returns the following matrix:
%
%      4    -1    -1     0
%     -1     4     0    -1
%     -1     0     4    -1
%      0    -1    -1     4
%
% See also: gallery, BuildLaplace1D, BuildLaplace2DStretchedMesh

      % N is the total number of points
      N=n^2;
      if exist('periodic') ~= 1
        periodic = false;
      elseif isa(periodic,'logical') == 0
        periodic = false;
      end

      fprintf('\n*****************************************************************************\n\n');
      fprintf('Building %d x %d scalar matrix\n',N,N);

      newMap = Map(N, 1);

      nsoffset = n*(n-1);
      ewoffset = n-1;
      MatrixData = spalloc(N,N,4*N);
      i=1; j=1; for k=1:N
        % Center point
        MatrixData(k,k) = 4;
        % North point
        if i ~= 1,  MatrixData(k,k-n) = -1;
        elseif periodic
          MatrixData(k,k+nsoffset) = -1;
        end
        % South point
        if i ~= n,  MatrixData(k,k+n) = -1;
        elseif periodic
          MatrixData(k,k-nsoffset) = -1;
        end
        % East point
        if j ~= 1,  MatrixData(k,k-1) = -1;
        elseif periodic
          MatrixData(k,k+ewoffset) = -1;
        end
        % West point
        if j ~= n,  MatrixData(k,k+1) = -1;
        elseif periodic
          MatrixData(k,k-ewoffset) = -1;
        end

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
      Amat = Operator(MatrixData,newMap,newMap,@MatlabApply, 'Laplace2D');

% debug test
%      i=42; M = BuildLaplace2D(i); M.MatrixData-gallery('poisson',i)
