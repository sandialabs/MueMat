%% A tutorial for new MueMat users
% In this example, we show how to use MueMat to solve the linear system $Ax=b$
% arising from the two-dimensional Laplace's equation, $\Delta u = f$, discretized
% on a square mesh with homogeneous Dirichlet boundary conditions.
%
% The linear system will be solved several times, each time with increasing levels
% of detail.

%% Preliminary remark
% MATLAB support two kinds of classes: handle classes and value
% classes. MueMat uses *handle classes*.
%
% Use of handle classes are natural if you are familiar with the
% concept of references and the oriented object support of other
% programming languages like C++ but please note that it contrasts
% with the rest of the Matlab language which is consistently
% pass-by-value.
%
% When you instantiate an handle class, the constructor returns a
% *reference* to the object created. When you assign this reference
% to another variable or pass it to a function, there is no copy of
% the original object and you still access to the initial data. In
% addition, a function that modifies the object does not need to
% return it.
%
% Example:
%
%   MgHierarchy = Hierarchy();     % Creation of an Hierarchy.
%   MgHierarchy2 = MgHierarchy;    % MgHierarchy2 is not a copy of MgHierarchy. It is just another name for accessing to the same data.
%   MgHierarchy.SetOutputLevel(1); % This method call also modifies MgHierarchy2.
%
% Most of MueMat classes implement a copy constructor and a Copy()
% method.
%
% Example:
%
%  MgHierarchy2 = Hierarchy(MgHierarchy); % MgHierarchy2 is a copy of MgHierarchy.
%  MgHierarchy3 = MgHierarchy.Copy();     % MgHierarchy3 is a copy of MgHierarchy.
%
% See also 'Comparing Handle and Value Classes' in the MATLAB
% documentation for more details about handle classes and value classes.

%% Example header
% We first initialize the random number generator to ensure reproducible results.
% The file |mue_include.m| defines some useful constants.
srand;
clear all;
mue_include;

%% Example 1: Minimimal information multigrid setup/solve:
% # Build the fine level matrix.
% # Allocate a smoothed aggregation data bucket and associate it with the fine level matrix.
% # Allocate the multigrid hierarchy and associate the smoothed aggregation data bucket created in the previous step with the finest level.
% # Populate the multigrid hierarchy using a smoothed aggregation factory.
% # Set the default smoothers.
% # Invoke the multilevel scheme either as solver or preconditioner.
%
% *Note*: At this point |MgHierarchy.Iterate()| does a fixed number of iterations.
% That is, there is no way to set a tolerance and iterate until the
% tolerance is met. It is recommended that |Iterate()| be used as a
% preconditioner.
%

Amat        = BuildLaplace2D(30);                            % Build the fine level matrix.
Finest      = Level();                                     % Allocate a smoothed aggregation data bucket.
Finest.Set('A', Amat);                                           % Associate the fine level matrix with the data bucket.
MgHierarchy = Hierarchy();                                         % Allocate the multigrid hierarchy
MgHierarchy.SetOutputLevel(1);                               % Verbose
MgHierarchy.SetLevel(Finest,1);                             % Associate the data bucket with the finest level.
MgHierarchy.FillHierarchy(SaPFactory());   % Populate the hierarchy using a smoothed aggregation factory.
MgHierarchy.SetSmoothers();                                  % Set the default smoothers.

n     = Amat.GetRowMap().NDOFs();
rhs   = rand(n,1);
guess = zeros(n,1);
CGtol = 100.*eps(2);
CGits = 99;

sol = MgHierarchy.Iterate(rhs, 9, guess, ALLZEROS);                         % Use multigrid as a solver.

[sol,flag,relres,iter,resvec]= pcg(Amat.GetMatrixData(),rhs,CGtol,CGits,...  % Use multigrid as a preconditioner.
                   @(rhs)MgHierarchy.Iterate(rhs,1, zeros(n,1), ALLZEROS));

%% Example 2: Understanding Iterate(rhs,Nits, guess, InitGuessStatus, Cycle, k)
% Optional arguments include
%
% * *|InitGuessStatus|*     |ALLZEROS| or |NOTALLZEROS| indicates whether the guess is all zeros or not. Some computations are skipped when |ALLZEROS|.  _Default_: |NOTALLZEROS|.
% * *|Cycle|*               1 ==> Vcycle, 2 ==> Wcycle.  _Default_: Vcycle.
% * *_k_*                 Iteration cycle starts on _k_ th level.  Normally, _k_ = 1 is the finest level.  _Default_: 1.

fprintf('\nSame solve as before but less efficient as initial guess is not\n');
fprintf('assumed to be the zero vector.\n');
sol = MgHierarchy.Iterate(rhs, 9, guess, NOTALLZEROS);
fprintf('\nWcycle solve\n');
sol = MgHierarchy.Iterate(rhs, 9, guess, ALLZEROS,2);
fprintf('\nNo output Vcycle solve\n');
sol = MgHierarchy.Iterate(rhs, 9, guess, ALLZEROS,1);

%% Example 3: Understanding SetSmoothers(Smfact, startLevel,numDesiredLevels, CrossFactory)
% All arguments are optional.
%
% *|Smfact|*
% Factory specifying smoothers to be built.
% Right now MueMat has only one smoother
% factory called |BlkSmootherFactory()| which
% itself takes arguments to specify things
% like Jacobi, GaussSeidel, BlkJacobi,
% BlkGaussSeidel, and domain decomposition.
% See |example/SmoothingTest.m| for a detailed
% description of |BlkSmootherFactory()|.
% _Default_:
% |BlkSmootherFactory('GaussSeidel', 1, 1)|,
% which corresponds to one pre and one post
% symmetric Gauss-Seidel iteration.
%
% *|startLevel|*          _Default_: 1.
%
% *|numDesiredLevels|*    _Default_: # of levels in hierarchy.
% |startLevel| and |numDesiredLevels| indicate
% which levels will have smoothers set.
% Note: By default |Smfact()| does not
% actually populate the last level with a
% smoother. A direct level solver is used
% on the coarsest level if no smoother is
% set. Thus, in reality |numDesiredLevels-1|
% smoothers are set with the default
% behavior of |Smfact|. This is a bit confusing
% but is consistent with |FillHierarchy()|.
%
% *|CrossFactory|*
% Handles any smoother request which
% involves coordination with a
% nonsmoother factory. This is passed to
% |Smfact|'s build method. An example
% of this might be using a block smoother
% where blocks are defined by aggregates
% generated from the prolongator factory.
% See |example/CrossFactorySpecTest.m| for
% a further discussion.
% _Default_: []
fprintf('\nV(3,3) cycle using Jacobi with omega= .7 !\n');
NewHierarchy = MgHierarchy.Copy(); % Copy the previously built hierarchy.
NewHierarchy.SetSmoothers(SmootherFactory(Smoother('Jacobi', 3, .7)));
sol = NewHierarchy.Iterate(rhs, 9, guess, ALLZEROS);
fprintf('\nControl of smoothers is described in example/SmoothingTest.m\n');
fprintf('Here is an example corresponding to a pre-forward Gauss-Seidel\n');
fprintf('and post-backward sweep Gauss-Seidel.\n');
NewHierarchy = MgHierarchy.Copy(); % Copy the previously built hierarchy.
PreSmoo     = Smoother('GaussSeidel', 1, 1);
PostSmoo    = Smoother(PreSmoo); % Copy
PreSmoo.SetForwardSweep(true);   PreSmoo.SetBackwardSweep(false);
PostSmoo.SetForwardSweep(false); PostSmoo.SetBackwardSweep(true);
SFact = SmootherFactory(PreSmoo,PostSmoo);
NewHierarchy.SetSmoothers(SFact);
sol = NewHierarchy.Iterate(rhs, 9, guess, ALLZEROS);
fprintf('\nHere we use the smoother on the coarsest level\n');
NewHierarchy = MgHierarchy.Copy(); % Copy the previously built hierarchy.
NewHierarchy.SetSmoothers(SFact);
NewHierarchy.SetCoarsestSolver(SFact);
sol = NewHierarchy.Iterate(rhs, 9, guess, ALLZEROS);
fprintf('\nFinally an example with different smoothers on different levels\n');
NewHierarchy = MgHierarchy.Copy(); % Copy the previously built hierarchy.
NewHierarchy.SetSmoothers(SmootherFactory(Smoother('GaussSeidel',1,1)),1,1);
NewHierarchy.SetSmoothers(SmootherFactory(Smoother('Jacobi',3, .7)),2,1);
sol = NewHierarchy.Iterate(rhs, 9, guess, ALLZEROS);

%% Example 4: Understanding FillHierarchy(PRfact, Acfact, startLevel, numDesiredLevels, CrossFactory)
% All but the first argument are optional.
%
% *|PRfact|*              Factory specifying transfer operators to be built.
%
% *|Acfact|*              Factory specifying coarse discretization
% matrices to be built.
% _Default_: RAPFactory() which just does
% Petrov-Galerkin projection.
%
% *|startLevel|*          _Default_: 1.
% numDesiredLevels    _Default_: # of levels in hierarchy.
% startLevel and numDesiredLevels indicate
% which levels will have operators built.
%
% *|CrossFactory|*        Handles any request which involves
% coordination between factories. It also
% handles all reuse/save requests (as these
% are sometimes cross factory).
% See example/CrossFactorySpecTest.m for
% a further discussion.
% _Default_: []


fprintf('Finally, we change a few default options corresponding to\n');
fprintf('a rectangular aggregation algorithm and lowering the \n');
fprintf('smoothed aggregation damping parameter\n');
CoalesceDropFact    = CoalesceDropFactory();
AggFact             = AggregationFactory();
AggFact.SetAlgorithm('rectangle');
AggPtsPerDim(1)     = 6; AggPtsPerDim(2) = 3;
AggFact.SetAggPtsPerDim(AggPtsPerDim);
AggFact.SetTargetSize(5);
PtentFact           = TentativePFactory(CoalesceDropFact,AggFact);
PFact               = SaPFactory(PtentFact);
PFact.SetDampingFactor(1./3);
PRFact              = GenericPRFactory(PFact);
MgHierarchy  = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PFact);
MgHierarchy.SetSmoothers();
sol = MgHierarchy.Iterate(rhs, 19, guess, ALLZEROS);
