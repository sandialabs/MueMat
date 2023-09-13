
srand;
clear all;
mue_include;

%% define problem
kappa = 0.01; % diffusivity constant within diffusion-convection eq (0.01 -> only mildly nonsymmetric)
sc = Scatra2D(20,kappa); % initialize doubleglazing problem with 20 x 20 elements
[Amat,rhs,problem] = sc.Build(); % generate problem (matrix, rhs, information about problem definition)
n     = Amat.GetRowMap().NDOFs();
guess = zeros(n,1);

%% prepare AMG
Finest      = Level();                          % Allocate a smoothed aggregation data bucket.
Finest.KeepAll(false);
Finest.Set('A', Amat);                                     % Associate the fine level matrix with the data bucket.

%% setup prolongation and restriction strategy
Pfact = SaPFactory();           % use smoothed aggregation (works, since problem only "mildly nonsymetric)
Rfact = GenericRFactory(Pfact);
PRfact = GenericPRFactory(Pfact,Rfact);

%% set smoother factory
Sfact = SmootherFactory(Smoother('Jacobi',3,0.66));

%% fill multigrid hierarchy
MgHierarchy = Hierarchy();                                   % Allocate the AMG hierarchy
MgHierarchy.SetOutputLevel(10);                         % Verbose
MgHierarchy.SetLevel(Finest,1);                       % Associate the data bucket with the finest level.
MgHierarchy.FillHierarchy(PRfact,RAPFactory(),1,3);
MgHierarchy.SetSmoothers(Sfact);                            % Set the default smoothers.


%% solve problem
sol = MgHierarchy.Iterate(rhs, 100, guess, ALLZEROS);                         % Use AMG as a solver.

% alternatively: use AMG as a precondtioner within GMRES
%[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,size(Amat.GetMatrixData(),1),1e-7,size(Amat.GetMatrixData(),1),@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

%% plot solution
sc.plot_vector(sol);
