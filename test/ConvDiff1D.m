%
% Simple 1D convection-diffusion (singular) problem with large diameter aggregates
% 
% Note: this differs from Penn. State slides in that it uses matlab's GMRES
%       with a W cycle.

NumIts = 4;
APdegree = 4; AggDiam = APdegree*2+1; Ncoarse = 99;  
Nfine = AggDiam*Ncoarse; 
h = 1/(Nfine+1);
%ii = input('Convection: 1=constant, 2=one jump, 3=smooth, 4=random \n');
ii = 3;
if (ii == 1) 
   c = input('value of c in u_xx + c u_x ');
   c = c*ones(Nfine,1);
elseif (ii == 2)
   clower = input('lower value of c in u_xx + c u_x ');
   cupper = input('upper value of c in u_xx + c u_x ');
   c = clower*ones(Nfine,1);
   c( floor(Nfine/2):Nfine ) = cupper;
elseif (ii == 3)
   c = .1*sin(2*pi*(1:Nfine)/Nfine)/h;
else
   c = rand(Nfine,1)/h;
end;
%
% build matrix
%
A = sparse(Nfine,Nfine);
for i=1:Nfine
    A(i,i) = 2/(h*h);
    if (i ~= Nfine) 
       A(i,i+1) = c(i)/h -1/(h*h);
    else
       A(i,1) = c(i)/h -1/(h*h);
    end;
    if (i ~= 1) 
       A(i,i-1) = -c(i)/h -1/(h*h);
    else
       A(i,Nfine) = -c(i)/h -1/(h*h);
    end;
end;



mue_include;

numDesiredLevels = 3;         % number of multigrid levels

AggFact  = AggregationFactory();
AggFact.SetAlgorithm('Uniform1D');
AggFact.SetAggPtsPerDim(AggDiam);
Patfact  = AP_PatternFactory();
Patfact.SetPatternType('AP');
Patfact.SetDegree(APdegree);
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact);
Pfact    = EminPFactory(Patfact, ConstraintFactory(), GMRESEminSolver(NumIts), Pinitfact); 
AltPfact = EminPFactory(Patfact, ConstraintFactory(), GMRESEminSolver(NumIts), Pinitfact);

Rfact    = GenericRFactory(AltPfact);
Acfact   = RAPFactory();
GSFactory= SmootherFactory(Smoother('GaussSeidel', 1, 1.)); 

%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.KeepAll(true);
Finest.Set('A', Operator(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
status = MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory, 1, numDesiredLevels);
nnn = size(A,1);
rhs = A*sin(2*pi*(1:nnn)/nnn)';

[sol,flag,relres,iter,resvec]=gmres(A,rhs,[],1e-11,100,@(v)MgHierarchy.Iterate(v,1,zeros(nnn,1),ALLZEROS,2),[],zeros(nnn,1));
fprintf('Nits = %d, ||r|| = %e \n', length(resvec)-1,norm(rhs-A*sol)/norm(rhs));

