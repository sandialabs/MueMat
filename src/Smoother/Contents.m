% SMOOTHER
%
% This directory contains smoother classes.
%
% Superclasses:
%
%   SmootherBase             - Base class for smoothers
%   SmootherPrototype        - Base class for smoother prototypes
%   SmootherFactoryBase      - Base class for smoother factories
%
% Smoother Factories:
%
%   SmootherFactory          - Generic Smoother Factory for generating the smoothers of the MG hierarchy
%   Hybrid2x2SmootherFactory - Specialized smoother factory for hybrid smoothers
%
% Available Smoothers:
%
%   Smoother                 - Jacobi, Gauss-Seidel smoother class (and more)
%   ChebySmoother            - Chebyshev smoother class
%   ILUSmoother              - ILU smoother class
%   Hybrid2x2Smoother        - Hybrid smoother: merges two basic smoothers for use on a 2x2 system
%   DirectSolveSmoother      - Direct solve smoother class
%
% Miscellaneous:
%
%   MakeUpRandomBlks         - Utility function for testing things
