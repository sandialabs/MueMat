%
% Test which runs a Simple-like smoother within a multigrid method applied
% to a 2x2 block system of the form
%             ( F  Bt )
%             ( B   0 )
%
%
% The smoother basically looks like
%
%[I      -GS(F) Bt]  [   GS(F)           0          ]  [  F         Bt ]
%[          I     ]  [ -iS B GS(F)       Cheby(Shat)]  [  B          C ]
%
%where GS(F) is running Gauss-Seidel (not symmetric GS) on F, and Cheby is
%running Chebyshev as a smoother. In the left results, Shat is
%implemented in a matrix-free way and corresponds to
%
%              Shat = C - B GS(F) Bt
%
%where an explicit diagonal is given to precondition the Chebyshev iteration
%and corresponds to  diag( C - B 1./rowsum(abs(F)) Bt ). In the right results,
%
%              Shat = C - B 1./(.3*rowsum(abs(F))) Bt ).
%
% Note: .3 scaling is just a HACK which sometimes works well.

function SIMPLESmooTest()

SetHomeDir
load([MUEMAT_ROOT_DIR '/data/SIMPLESmoother.mat']);

nv   = size(F,1);
np   = size(B,1);
pres = (3:3:nv+np)';
vel  = ones(nv+np,1); vel(pres) = 0; vel = find(vel);
rhs  = rhs([vel; pres]);
A    = [F Bt; B C];
Lap  = -spones(C);  Lap = Lap-spdiags(sum(Lap,2),0,np,np);

mue_include;

AggFact   = AggregationFactory();
FCFact    = FCSplittingFactory(AggFact);
Pinitfact = TentativePFactory(CoalesceDropFactory,AggFact);
PatFact   = AP_PatternFactory([], Pinitfact);

EFact = EminPFactory(PatFact, ConstraintFactory(PatFact, FCFact), [], Pinitfact);
EFact.SetEminSolverFunction(GMRESEminSolver(3));
Pfact   = PFactoryForSimple(EFact,SaPFactory());
Rfact   = TransPFactory();
Acfact  = RAPXfemFactory(FCFact);

Finest = Level();   
Finest.KeepAll(false);       % SimpleSmoother needs special variables, supports no SetNeeds -> Keep them
Finest.Keep('N11'); Finest.Keep('A11');
Finest.Keep('N22'); Finest.Keep('A22');
Finest.Keep('NullSpace');
Finest.Set('A',   Operator(A,   Map(nv+np,1),Map(nv+np,1),@MatlabApply));
Finest.Set('A11', Operator(F,   Map(nv/2,2), Map(nv/2,2), @MatlabApply));
Finest.Set('A22', Operator(Lap, Map(np,1),   Map(np,1),   @MatlabApply));
Finest.Set('N11', size(Finest.Get('A11'),1));
Finest.Set('N22', size(Finest.Get('A22'),1));

temp = zeros(nv,2); temp(1:2:end,1)=1; temp(2:2:end,2)=1;
Null.top = temp;
Null.bot = ones(np,1);
Finest.Set('NullSpace', Null);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, NLevels);

Direct = DirectSolveSmoother(); 
GS     = Smoother('GaussSeidel', NFSweeps, 1.);
GS.SetForwardSweep(true);
GS.SetBackwardSweep(false);
Cheb   = ChebySmoother(NSSweeps, 1/9.,'mydiag');

Simp = SimpleSmoother(GS,GS,Cheb,1);
if UseDinv, Simp.SetUseDinv(true); end; 
SFact = SmootherFactory(Simp,[]);  % No post smoothing in this example
MgHierarchy.SetSmoothers(SFact);
MgHierarchy.SetCoarsestSolver(SFact);

params(1) = 1e-8; params(2) = 200; x0   = zeros(nv+np,1);
[sol,flag,relres,iter,resvec]=gmres(A,rhs,[],1e-8,200,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(nv+np,1),ALLZEROS),[],x0);
format long e; fprintf('# its = %d\n',length(resvec)-1); resvec,
end
