% function [Amat, NullSpace] = BuildElasticity2D(n, E, nu,[stretch])
% Build a matrix (and near-nullspace) for 2D linear elasticity.
%
% This uses a tri-linear elements on a regular mesh (allows for
% mesh stretching, if requested)
%
% Parameters:
% 'n'  - Number of grid points in each direction.
% 'E'  - Elastic modulus.
% 'nu' - Poisson's ratio.
% Optional Parameters:
% 'stretch' - A 2-vector containing the mesh streching in each
%             dimension ([1,1] = isotropic mesh).
%
% Example:
% [M,nullspace]=BuildElasticity2D(2 ,1e5, 0.3);
%
% References:
% [1] R.D. Cook, D.S. Malkus and M.E. Plesha
%     "Concepts and Applications of Finite Element Analysis, 3rd ed.
%
% See also: BuildElasticity3D BuildElasticity

function [Amat, NullSpace,NODES,ELEMENTS] = BuildElasticity2D(n, varargin)
    if ~(nargin>1 && strcmp(varargin{nargin-1},'old'))
      [Amat, NullSpace, NODES,ELEMENTS] = BuildElasticity([n+1,n+1], varargin{:});
    else
      % to get the old Elasticity2D test case, call
      % BuildElasticity2D(n, E, nu, 'old')
      fprintf('Elasticity2D: Running Old\n');
      if(nargin>2), E=varargin{1}; else E=1e9;end
      if(nargin>3), nu=varargin{2}; else nu=.25;end
      [Amat, NullSpace, NODES] = BuildOldElasticity2D(n, E, nu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another version of BuildElasticity2D follows
% - This version return the same matrix as the 2D elasticity test case of pyAMG
% - But no mesh streching ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build the matrix and the near nullspace matrix of a 2D Linear Elasticity problem
%
% This example is a matlab implementation of an example of PyAMG
% (http://code.google.com/p/pyamg/)
%
% It is a finite difference discretization (Q1 finite elements) of
% a linear elasticity problem on a grid.
%
% Parameters:
% 'n'  is the number of grid point on each direction.
% 'E'  is the Young's modulus.
% 'nu' is the Poisson's ratio.
%
% The function return a scalar sparse matrix (Amat) and the
% corresponding nullspace (as a matlab matrix).
%
% Example:
% [M,nullspace]=BuildElasticity2D(2 ,1e5, 0.3); M.MatrixData
%
% References:
% [1] N. Bell, L. Olson and J. Schroder
%     PyAMG implementation of the same problem:
%     pyamg/gallery/elasticity.py, function q12d
% [2] J. Alberty, C. Carstensen, S. A. Funken, and R. KloseDOI
%     "Matlab implementation of the finite element method in
%     elasticity" - Computing, Volume 69, Issue 3 (Nov. 2002)
%     Pages: 239 - 263
%     http://www.math.hu-berlin.de/~cc/
%
% See also: BuildMatrixElasticity2D BuildNullSpaceElasticity2D
%
function [Amat, NullSpace,NODES] = BuildOldElasticity2D(n, E, nu)
      dirichlet_boundary = true;
      % dirichlet_boundary = false;
      if dirichlet_boundary
        n=n+1;
      end

      MatrixData = BuildMatrixElasticity2D(n,E,nu);
      [NullSpace,NODES] = BuildNullSpaceElasticity2D(n);

      % Add Dirichlet boundary condition
      if dirichlet_boundary
        P=AddBoundaryCondition(n);
        MatrixData = P.' * MatrixData * P;
        NullSpace = P.' * NullSpace; % simplify ?
        clear P
      end

      % Fill-in Amat members
      X=(n-1)^2; % X=(n+1)^2; if dirichlet_boundary=false
      newMap = Map(2*X, 1);
      %Amat.RowMap     = newMap;
      %Amat.ColMap     = newMap;
      %Amat.Apply      = @MatlabApply;
      %Amat.MatrixData = MatrixData;
      Amat = Operator(MatrixData,newMap,newMap,@MatlabApply, 'Elasticity2D');

      NullSpace=normalize(NullSpace);


% Build the global matrix of a 2D linear elasticity problem.
% The function return a scalar sparse matrix (matlab matrix).
% This function is a matlab implementation of an example of PyAMG
%
% The function return a scalar sparse matrix (matlab format).
%
% See BuildElasticity2D for more details.
%
function M = BuildMatrixElasticity2D(n, E, nu)
      N = n^2;

      lame = E * nu / ((1 + nu) * (1 - 2*nu)); % Lame's first parameter
      mu   = E / (2 + 2*nu);                   % Shear modulus

      % Compute local stiffness matrix (size=8x8)
      K = BuildLocalMatrixElasticity2D(lame, mu);

      % Build the global stiffness matrix (format I,J,V)
      % a) Mesh grid (LL)
      %    if n=2, LL=[0 4 8; 1 5 9; 2 6 10]
      nodes = reshape(meshgrid(0:(n+1)^2-1,1),n+1,n+1);
      LL =  nodes(1:n,1:n);

      % b) Index of matrix coefficients
      %    Note: index begin with 0 (we add +1 later)
      I = 2*LL; % I = 2*LL+1;
      I = reshape(I,1,[]);
      I = reshape(repmat(I.', 1, 8*8).',1,[]);
      J = I; % here, for n=2, I=J=[0...0 2...2 8...8 10...10 12...12 etc.]

      mtx = repmat([0, 1, 2, 3, 2*n + 4, 2*n + 5, 2*n + 2, 2*n + 3],8,1); %tmp
      I = I + reshape(repmat(mtx.', 1,N),1,[]);
      J = J + reshape(repmat(mtx  , 1,N),1,[]);
      clear mtx;

      % c) Set values of coefficients
      %    V = [ K K K ... K ]
      V  = repmat(reshape(K,1,[]), 1, N);

% debug
%      disp(I); disp(J); disp(V);
%      save('test', 'I', 'J', 'V', '-ascii');

      % Convert (I,J,V) to a matlab sparse matrix
      % (index begin with 1)
      M = sparse(I+1,J+1,V);
      clear I J V;


% Build the local matrix of a 2D Linear Elasticity problem
% on a square element (Q1)
% This function is a matlab implementation of an example of PyAMG
%
% Parameters:
% 'lame' is the Lame's first parameter
% 'mu'   is the Shear modulus
%
% The function return a scalar sparse matrix of size 8x8
% (matlab format).
%
% Vertices should be listed in counter-clockwise order:
%
%         [3]----[2]
%          |      |
%          |      |
%         [0]----[1]
%
% Degrees of freedom are enumerated as follows:
%
%    [x=6,y=7]--[x=4,y=5]
%        |          |
%        |          |
%    [x=0,y=1]--[x=2,y=3]
%
function K = BuildLocalMatrixElasticity2D(lame, mu)

      M = lame + 2*mu; % P-wave modulus

      R_11 = ([2 -2 -1  1;
               -2  2  1 -1;
               -1  1  2 -2;
               1 -1 -2  2]) / 6.0;

      R_12 = ([1  1 -1 -1;
               -1 -1  1  1;
               -1 -1  1  1;
               1  1 -1 -1]) / 4.0;

      R_22 = ([2  1 -1 -2;
               1  2 -2 -1;
               -1 -2  2  1;
               -2 -1  1  2]) / 6.0;

      K=zeros(8,8);

      E = [M 0; 0 mu];
      K([1 3 5 7],[1 3 5 7]) = E(1,1) * R_11 + E(1,2) * R_12 + E(2,1) * R_12' + E(2,2) * R_22;

      E = [mu 0; 0 M];
      K([2 4 6 8],[2 4 6 8]) = E(1,1) * R_11 + E(1,2) * R_12 + E(2,1) * R_12.' + E(2,2) * R_22;

      E = [0 mu; lame 0];
      K([2 4 6 8],[1 3 5 7]) = E(1,1) * R_11 + E(1,2) * R_12 + E(2,1) * R_12.' + E(2,2) * R_22;

      K([1 3 5 7],[2 4 6 8]) = K([2 4 6 8],[1 3 5 7]).';


% Build a matrix 'filter' to add Dirichlet boundary condition to
% the 2D elasticity problem.
%
function P = AddBoundaryCondition(n)
        % 'mask' is used to find boundary of the mesh (n+1 x n+1):
        % for n=3, mask =
        %      0     0     0     0
        %      0     1     1     0
        %      0     1     1     0
        %      0     0     0     0
        mask=zeros(n+1,n+1);
        mask(2:n, 2:n) = 1;

        % matrix format: i,j,v
        % compute (i,j) indexes corresponding to boundaries on the
        % stiffness matrix
        i = 2*find(mask).'-1;
        j = 2*meshgrid(1:(n-1)^2,1)-1;

        i = [i i+1]; % There are 2 coefficients per mesh point.
        j = [j j+1];

% debug
%       disp(i); disp(j);

        v = ones(1, (n-1)^2);
        v = [v v];

        P=sparse(i,j,v, 2*(n+1)^2, 2*(n-1)^2);

% debug
%      full(P)


% Build the nullspace matrix for the 2D Linear Elasticity problem
% on a regular grid.
% There are 3 near nullspace modes.
%
% See also: BuildElasticity2D build_elastic_rbm
%
function [NullSpace,NODES] = BuildNullSpaceElasticity2D(n)
      x=n+1;

      % 2D mesh
      NODES(:,1) = reshape(meshgrid(1:x,1:x).', 1, []).';  % = [ 1234 1234 1234 1234 ]
      NODES(:,2) = reshape(meshgrid(1:x,1:x),   1, []);    % = [ 1111 2222 3333 4444 ]

      % nullspace is the rigid-body mode matrix
      NullSpace = build_elastic_rbm(NODES);

% Given a set of nodes in 1, 2 or 3 dimensions construct the
% rigid-body modes for the (elastic) structure
%
% by: Chris Siefert
% Last Updated: 11/02/2006
%
% See also: BuildNullSpaceElasticity2D
%
function RBM=build_elastic_rbm(NODES)
      % Constants
      [N,dim]=size(NODES);
      NDOF=[1,3,6];

      % Allocs
      RBM=zeros(dim*N,NDOF(dim));

      % Translational DOFs / INDICES
      for I=1:dim,
        IDX{I}=I:dim:dim*N;
        RBM(IDX{I},I)=ones(N,1);
      end

      % Recenter nodes
      CTR=sum(NODES) / N;
      CNODES = NODES - repmat(CTR,N,1);

      % Rotational DOF:  Thanks to Farhat, Pierson and Lesoinne 2000,
      % equation (52), you know, once I've managed to make sense out of
      % the blasted thing.
      if(dim>=2)
        % Rotate in X-Y Plane (around Z axis): [-y ;x];
        RBM(IDX{1},dim+1)=-CNODES(:,2);
        RBM(IDX{2},dim+1)= CNODES(:,1);
      end
      if(dim==3)
        % Rotate in Y-Z Plane (around X axis): [-z;y]
        RBM(IDX{2},dim+2)=-CNODES(:,3);
        RBM(IDX{3},dim+2)= CNODES(:,2);

        % Rotate in X-Z Plane (around Y axis): [z ;-x]
        RBM(IDX{1},dim+3)= CNODES(:,3);
        RBM(IDX{3},dim+3)=-CNODES(:,1);
      end

      RBM=normalize(RBM);


% Happy normalize function (normalizes columns)
function [AN,NF]=normalize(A)
SZ=size(A,2);
NF=zeros(SZ,1);
AN=0*A;
for I=1:SZ,
  NF(I)=norm(A(:,I));
  if(NF(I) < 1e-11)
    NF(I)=1;AN(:,I)=A(:,I);
  else
   AN(:,I)=A(:,I)/NF(I);
  end
end
