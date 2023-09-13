%% A first tutorial on the use of the facade classes
% This small tutorial demonstrates the usage of a small MueMat facade class
% for elasticity problems
%% Idea
% The idea of facade classes is to provide an easy-to-use interface to
% MueMat that helps new MueMat users with their first steps in MueMat.

%% Example: elasticity problem
% First we have to load the problem matrix and the corresponding nullspace
% information. Note, that MueMat has only routines for creating default
% nullspaces with constant vectors. For elasticity problems the user itself
% has to provide nullspace information (rigid body modes) for his problem.

%%
% Our example is a simple 2d elasticity problem. We load the data file with
% the matrix and the nullspace information. Since it's only 2d we have a
% three dimensional null space (translation in x and y direction and the
% rotational degrees of freedom).

clear;
mue_include;

SetHomeDir
mytests = { [MUEMAT_ROOT_DIR '/data/TutorialFacade.mat']};
load(mytests{1});

%%
% The data file only contains a MATLAB matrix. For MueMat we need a MueMat
% operator object. We transform the MATLAB matrix to the MueMat operator
% object with 2x2 block maps (2 dimensional problem!) using

A = Operator(A,2,2);

%%
% Since it's an elasticity problem, we use the facade class for elasticity
% problems that is called |Elasticity_AMG|. Internally it uses symmetric
% SA-AMG (smoothed aggregation) transfer operators. The facade class has
% only on static member function |Setup|, that we can call for creating the
% multigrid hierarchy for a given problem matrix _A_ and the corresponding
% nullspace _nsp_
MgHierarchy = Elasticity_AMG.Setup(A,nsp);

%%
% Finally we use the multigrid hierarchy as a preconditioner for the
% standard MATLAB GMRES routine (of course you also can use PCG or other
% iterative solvers, or the Iterate routine of MueMat itself, if you want
% to use a plain Multigrid solver).
GMREStol = 1e-10;
GMRESits = 50;
[x,flag,relres,iter,resvec]= gmres(A.GetMatrixData(),b,A.GetRowMap().NDOFs(),GMREStol,GMRESits,...
        @(b)MgHierarchy.Iterate(b,1,zeros(length(b),1),ALLZEROS,1));
fprintf('flag = %d; res = %6.3g after %d iterations \n',flag,resvec(end),length(resvec)-1);

%% conclusion
% This short tutorial shows the usage of a small facade class for
% generating multigrid hierarchies in MueMat for elasticity problems. With
% only one row of code you're able to create a multigrid
% solver/preconditioner that you can apply to your problems.