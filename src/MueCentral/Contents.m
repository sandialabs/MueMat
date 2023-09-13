% MUECENTRAL
%
% Superclasses:
%
%   CopiableHandle      - Superclass to add copy capability to Matlab handle classes
%   FactoryBase         - Base class for factories
%
% Multigrid algorithm:
%
%   Hierarchy           - Provides methods to build a multigrid hierarchy and apply multigrid cycles
%   RAPFactory          - Factory which builds coarse discretizations via R*A*P
%
% Storing data of MG levels:
%
%   Level - Abstract class for data associated with a single level of a multigrid hierarchy
%   Level        - Data associated with one level of a smoothed aggregation method
%
% Interfactories communication:
%
%   CrossFactory        - Handles parameter specifications that affect several factories
%
% Miscellaneous:
%
%   mue_include         - MueMat header file
%
%  Hierarchy - Hierarchy Provides methods to build an multigrid hierarchy and apply multigrid cycles
%  SingleLevelFactoryBase - SINGLELEVELFACTORYBASE Abstract factory for building objects that require only a single level.
%  VerboseObject - Base class to control output. All other classes derive from this one.
