%% Demonstration of different ways to construct smoothers.
%
% This example builds a 1D Laplacian and solves via a 2-level multigrid hierarchy.
% The example shows different ways smoothers are constructed. In particular,
%
%  sym Gauss-Seidel       | SFact= BlkSmootherFactory('GaussSeidel',   Nits,Omega);
%  -----------------------------------------------------------------------------
%  sym block Gauss-Seidel | SFact= BlkSmootherFactory('BlkGaussSeidel',Nits,Omega);
%  -----------------------------------------------------------------------------
%  Jacobi                 | SFact= BlkSmootherFactory('Jacobi',        Nits,Omega);
%  -----------------------------------------------------------------------------
%  block Jacobi           | SFact= BlkSmootherFactory('BlkJacobi',     Nits,Omega);
%  -----------------------------------------------------------------------------
%  additive domain        | SFact= BlkSmootherFactory('Jacobi',Nits,Omega,...
%  decomposition          |                        'Random NonOverlapping');
%  -----------------------------------------------------------------------------
%  multiplicative domain  | SFact= BlkSmootherFactory('GaussSeidel',Nits,Omega,...
%  decomposition          |                        'Random NonOverlapping');
%  -----------------------------------------------------------------------------
%
% *Remarks*:
%
% # Omega == 1 runs the fast version of Gauss-Seidel. No fast 
% version of Jacobi is currently available (for doing Dinv). 
% # #pre iterations = # post iterations = |Nits|
% # domains are silly (random) for domain decomposition. Overlapping
% is also available.  The example is meant to demonstrate the 
% interface and hopefully more realistic things will be created.
%
%
% Fine grain control is available. For example,
%
%       SFact = SmootherFactory(Smoother('GaussSeidel', Nits, Omega));
%       SFact.SetIts(Npre,Npost);
%       SFact.SetOmega(OmegaPre, OmegaPost);
%       SFact.ForwardSweeps(TrueOrFalse,TrueOrFalse);
%       SFact.BackwardSweeps(TrueOrFalse,TrueOrFalse);
%
% *Remarks*:
%
% # |SetIts()| over-rides |Nits| in constructor and allows the 
% number of pre- and post- iterations to differ.
% # |SetOmega()| over-rides |Omega| in constructor and allows the 
% damping factor for pre- and post- iterations to differ.
% # |SetIts()| and |SetOmega()| do similar things for other smoothers.
% # |ForwardSweeps()| and |BackwardSweeps()| decide whether forward
% and/or backward sweeps are used within Gauss-Seidel iterations.
% They do similar things for block Gauss-Seidel and 
% multiplicative domain decomposition.
%
%% Advanced Concepts
%
% # The last example is for setting domain decomposition blocks. The basic
% idea in factories is that no real problem specific data is given when 
% factories are constructed. This means that only strategies for defining
% blocks can be given during the factory construction phase.  The standard
% procedure for setting smoothers is to create a smoother factory as
% illustrated above and pass this to |MgHierarchy.SetSmoothers()| which 
% automatically invokes the smoother factory's build method to create
% smoothers on all levels.  If we happen to actually have block data for
% a particular level, we must directly call the factory's build method 
% so that we can pass this data in. See last example for more details.
% 
% # If the number of pre (or post) iterations is equal to zero, then
% a pre (or post) smoothing object is not created when the factory's
% build method is invoked. This means that we can build completely 
% independent pre and post smoothing objects by invoking
% |MgHierarchy.SetSmoothers()| twice with two different smoother
% factories. The first time setting PostIts to 0 and the second time
% setting |PreIts| to 0. 
%

%%
% *Problem setup*
clear;
srand;
mue_include
InitGuessStatus = NOTALLZEROS;

Amat = BuildLaplace1DBlk(1,-1,100);

rowmap = Amat.GetRowMap();
rhs  = rand(rowmap.NDOFs(),1);
sol  = zeros(rowmap.NDOFs(),1); InitGuessStatus = ALLZEROS;

%%
% *Factories to define AMG*
%
Pfact             = SaPFactory();
Rfact             = TransPFactory();
PRfact            = GenericPRFactory(Pfact,Rfact);
Acfact            = RAPFactory();
%%
% *Construct and populate finest level with user information*
%
Finest = Level();
Finest.KeepAll(false);
Finest.Set('A', Amat);
%Finest.Set('NullSpace', ones(rowmap.NDOFs(),1));
Finest.Set('NullSpace', BuildNullSpace(Amat));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 2);

%%
% *Different ways to set and exercise smoothers*
%
sweeps = 1;               % number of smoothing sweeps (pre or post)
startLevel = 1;           % start indexing at level 1
numDesiredLevels = 2;     % max number of levels in hierarchy
iterations = 9;           % multigrid iterations
smoothCoarsest = false;   % if true, run smoother on coarsest problem
                          % (must be true if a one-level method and you
                          %  don't want a direct solve)

%%
% *Symmetric Gauss-Seidel*
fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'point'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Symmetric Gauss-Seidel*
fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99 + smoother on coarsest problem\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', 3, .99, 'point'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (~smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, 9, sol, InitGuessStatus);

%%
% *Symmetric Block Gauss-Seidel, pre and post smoothing*
fprintf('Running symmetric Block Gauss-Seidel V(1,1) with w = .99\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'default'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Jacobi*
fprintf('Running Jacobi V(1,1) with w = .7\n');
SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7, 'point'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Block Jacobi*
fprintf('Running Block Jacobi V(1,1) with w = .7\n');
SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7, 'default'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Symmetric Gauss-Seidel, pre smoothing only*
fprintf('Running symmetric Gauss-Seidel V(1,0) with w = .99\n');
SFact    = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'point'), []);
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Symmetric Gauss-Seidel, post smoothing only*
fprintf('Running symmetric Gauss-Seidel V(0,1) with w = .99\n');
SFact    = SmootherFactory([],Smoother('GaussSeidel', sweeps, .99, 'point'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Symmetric Block Gauss-Seidel, pre and post smoothing but with different damping factors*
fprintf('Running pre-forward sweep only with omega = .5 and post-backward sweep only Gauss-Seidel V(1,1) with w = .99\n');
PreSmoo  = Smoother('GaussSeidel', sweeps);
PreSmoo.SetOmega(.5);
PreSmoo.SetForwardSweep(true);
PreSmoo.SetBackwardSweep(false);
PostSmoo  = Smoother('GaussSeidel', sweeps);
PostSmoo.SetOmega(.99); 
PostSmoo.SetForwardSweep(false); 
PostSmoo.SetBackwardSweep(true);
SFact    = SmootherFactory(PreSmoo,PostSmoo);
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Overlapping additive Schwarz, pre and post smoothing, random domains*
fprintf('Running overlapping additive domain decomp. V(1,1) where domains are chosen randomly\n');
SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7,'Random NonOverlapping'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Overlapping multiplicative Schwarz, pre and post smoothing, random domains*
fprintf('Running symmetric overlapping multiplicative domain decomp. V(1,1) where domains are chosen randomly\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99,'Random Overlapping'));
NewMg    = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *Alternative setup method for overlapping multiplicative Schwarz, pre and post smoothing, random domains*
fprintf('Running symmetric overlapping multiplicative domain decomp. V(1,1) where domains are chosen randomly but ...\n');
fprintf('setting the blocks in a more direct but low-level fashion\n');
NewMg = MgHierarchy.Copy();
Finest      = NewMg.GetLevel(1);
FineA       = Finest.Get('A');

%%
% *Switching views of matrix and smoothing using different diagonal*
FineA.CreateView('ovblock', FineA.GetRowMap(), FineA.GetColMap(), FineA.GetApply());
% create a new block diagonal
Collection = MakeUpRandomBlks(FineA,'Overlapping');
BlkDiag = FineA.GetDiagonal(Collection, 'ovblock');
if isempty(BlkDiag.GetApplyInverse())
  FactorBlkDiag(BlkDiag);
end

SFactLvl1 = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'ovblock'));
SFact     = SmootherFactory(Smoother('GaussSeidel', sweeps, .99));

NewMg.SetSmoothers(SFactLvl1, startLevel,   startLevel);
%NewMg.SetSmoothers(SFact, startLevel+1, numDesiredLevels-1); %TODO: useless: we need a 3 level problem here

if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;

newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);