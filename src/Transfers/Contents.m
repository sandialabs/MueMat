% TRANSFERTS
%
% This directory contains classes to build transfert operators
% between grids of the MG hierarchy.
%
% Superclasses:
%
%   TwoLevelFactoryBase - Abstract factory for building P, R, and A on coarse levels
%   PFactory            - Abstract interface factory for prolongator (Build)
%   RFactory            - Abstract interface factory for restrictor (Build)
%   PRFactory           - Abstract interface for prolongation and restriction transfer operator (Build)
%
% Standard transfer operators
%
%   TentativePFactory   - factory class for tentative prolongation operator (implements PFactory)
%   TentativeRFactory   - factory class for tentative restriction operator (implements RFactory)
%   TransPFactory       - Factory which builds restrictors by transposing prolongators
%   GenericPRFactory    - default implementation of PRFactory. puts PFactory object and RFactory object together ot PRFactory object.
%
% Smooth-aggregation (SA):
%
%   SaPFactory          - Build a prolongator via the Smoothed Aggregation algorithm (implements PFactory)
%
% Petrov-Galerkin smoothed aggregation (PG):
%
%   PgPRFactory         - PgPRFactory Build a prolongator and restrictor via the PG-AMG algorithm
%   PgPRFactory2        - same as PgPRFactory class. small optimization within code
%
% Energy minimization (Emin):
%
%   EminPFactory        - This factory creates a prolongator that is the solution of a constrained
%   CoarseNSFactory     - Factory for generating a coarse grid representation of the nullspace
%   EminPFactory2       - Build a prolongator via a constrained energy minimization algorithm.
%
% Energy minimization (MinDescent):
%
%   MinDescentPRFactory - This factory creates a prolongation operator that is improved by a constrained minimization method
%
% Sparsity Pattern factories:
%
%   PatternFactory           - abstract bactory class for patterns
%   AP_PatternFactory        - concrete implementation of PatternFactory for AP derived patterns
%   AffInvAfc_PatternFactory - concrete implementation of PatternFactory for AffInvAfc derived patterns
%
% Sparsity Pattern (Filtering):
%
%   PatternFilter            - abstract class for sparsity pattern filters (for prolongator and restrictor sparsity patterns)
%   ConstEntriesPerRowFilter - concrete implementation of PatternFilter (const number of nonzeros per row in prolongator)
%   Thresholding             - concrete implementation of PatternFilter (just thresholding)
%
% Miscellaneous:
%
%   GenericRFactory     - This restriction factory wraps a prolongation factory (deprecated?)

%  SaCoarseNSFactory - A factory for generating a coarse grid representation of the

%  PatternFactory_old - PATTERNFACTORY Factory to create prolongator sparsity patterns
%  AvgNNZperRowFilter - AvgNNZperRowFilter
%  ConstNNZperRowFilter - copy constructor
%  InlineNNZperRowFilter - InlineNNZperRowFilter constructor
%  PatternFilterFactory - PATTERNFILTER base class for a transfer operator pattern filter
%  RelativeThresholding - RELATIVETRESHOLDING implementation of a relative thresholding pattern filter
