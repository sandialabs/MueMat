% function [Amat, NullSpace] = BuildElasticity3D(n, E, nu,[stretch])
% Build a matrix (and near-nullspace) for 3D linear elasticity.
%
% This uses a tri-linear elements on a regular mesh (allows for
% mesh stretching, if requested)
%
% Parameters:
% 'n'  - Number of grid points in each direction.
% 'E'  - Elastic modulus.
% 'nu' - Poisson's ratio.
% Optional Parameters:
% 'stretch' - A 3-vector containing the mesh streching in each
%             dimension ([1,1,1] = isotropic mesh).
%
% Example:
% [M,nullspace]=BuildElasticity3D(2 ,1e5, 0.3);
%
% References:
% [1] R.D. Cook, D.S. Malkus and M.E. Plesha
%     "Concepts and Applications of Finite Element Analysis, 3rd ed.
%
% See also: BuildElasticity2D BuildElasticity

function [Amat, NullSpace,NODES] = BuildElasticity3D(n, varargin)
    [Amat, NullSpace, NODES]  = BuildElasticity([n+1,n+1,n+1], varargin{:});
%   [Amat, NullSpace, NODES]  = BuildElasticity([n,  n,  n], varargin{:});
end
