function GeoDDTest(doomversion,stretch,smoonxdomain,smoonydomain,nxdomain,nydomain,numDesiredLevels)
%  Should really fix things so that direct solver on coarsest grid is replaced
%  with something that projects out null space before inverting.
%
% Jeremie,
%    I wanted to check in another slightly strange example that I have been 
% working with. I want to make sure that I am doing things right with 
% respect to requests/gets and I wanted to point out a couple of unique
% aspects to this problem. 
%    1) In my first version, I made GeometricAggFact() to produce aggregates.
%       This factory doesn't really need anything dropped/coalesced, so
%       I changed EminPFactory so that if an empty factory was given, it
%       just skips the coalesce.build() phase. 
%    2) I did AggFact.SetGraphName('A'); because I just wanted the matrix A
%       used as the graph. I'm not sure if this is what I am supposed to do.
%       I'm not sure if something bad could happen with release.
%    3) I actually don't need an aggregation factory ... just a pattern
%       factory ... so I tried removing it. I decided not to continue with
%       this as it is tricky, and I'm hopeful that the return of BabyEmin
%       will make this easier.
%    4) I also wanted to modify the coarsest grid direct solve so that
%       it projects out the null space of the singular matrix before
%       and after the solve. I just put some code into Hierarchy.m that is
%       now commented out. I wasn't sure the best way to do this, so I
%       thought we could talk over the phone.

% This example defines a set of regular domains and cross-points on 
% an unstructured mesh and effectively does a 2-level domain decomposition
% method where the coarse grid unknowns are the cross-points and the
% sparsity pattern of the grid transfers is geometrically defined.

SetHomeDir

doom=1;

if exist('numDesiredLevels') ~= 1
  numDesiredLevels     = 2;
end
if numDesiredLevels == 1
  fprintf('\n** TESTING SMOOTHER ONLY **\n\n');
end

%stretch=100;

%doomversion=1;
%doomversion=2;
%doomversion=3;
%doomversion=4;
%doomversion=5;
%doomversion=6;
%doomversion=7;

%dataFile = 'newadapt';
%dataFile = 'newadapt8';
%dataFile = 'newadapt9';
%dataFile = 'newadapt10';

showmesh = false;

%smoonxdomain = 15; smoonydomain = 15;

%nxdomain = 4; nydomain = 4;
%nxdomain = 6; nydomain = 6;
%nxdomain = 8; nydomain = 8;
%nxdomain = 10; nydomain = 10;
%nxdomain = 12; nydomain = 12;
%nxdomain = 16; nydomain = 16;
%nxdomain = 16; nydomain = 16;
%nxdomain = 20; nydomain = 20;
%nxdomain = 26; nydomain = 26;
%nxdomain = 30; nydomain = 300;
%nxdomain = 50; nydomain = 50;
%nxdomain = 76; nydomain = 76;
%nxdomain = 100; nydomain = 100;

vizMesh=false;

overlap=1;
linespec='k-.';
maxIts = 100; 
stopTol = 1e-10;

if doom
  loadDataStr = ['[A,mesh] = tester(' num2str(doomversion) ',' num2str(stretch) ',showmesh);'];
else
  loadDataStr = ['load ' MUEMAT_ROOT_DIR '/data/' dataFile];
end

h = 2;
%h = 3;
initialH = 1/h; % width of coarsest domain
%
% setup the smoother subdomains
%
smooDomainList.Ndomains = smoonxdomain*smoonydomain;
smooDomainList.dim = 2;
smooDomainList.LowerLeftCorner = zeros(smooDomainList.Ndomains,smooDomainList.dim);
smooDomainList.UpperRightCorner= zeros(smooDomainList.Ndomains,smooDomainList.dim);
%smooDomainList.Bcs =  char(double('i')*ones(smooDomainList.Ndomains,smooDomainList.dim*2));
smoodx = 1/smoonxdomain;
smoody = 1/smoonydomain;

for i= 1:smoonxdomain
   for j= 1:smoonydomain
      smooDomainList.LowerLeftCorner((j-1)*smoonxdomain+i,1) = (i-1)*smoodx;
      smooDomainList.LowerLeftCorner((j-1)*smoonxdomain+i,2) = (j-1)*smoody;
      smooDomainList.UpperRightCorner((j-1)*smoonxdomain+i,1)= (i  )*smoodx;
      smooDomainList.UpperRightCorner((j-1)*smoonxdomain+i,2)= (j  )*smoody;
   end
end

for j= 1:smoonydomain
      smooDomainList.Bcs((j-1)*smoonxdomain+1,1) = 'b';
      smooDomainList.UpperRightCorner(j*smoonxdomain,1)= 1.000001;
      smooDomainList.Bcs(j*smoonxdomain,2) = 'b';
end
for i= 1:smoonxdomain
   smooDomainList.Bcs(i,3) = 'b';
   smooDomainList.UpperRightCorner((smoonydomain-1)*smoonxdomain+i,2)= 1.000001;
   smooDomainList.Bcs((smoonydomain-1)*smoonxdomain+i,4) = 'b';
end

dx = 1/nxdomain;
dy = 1/nydomain;
%
%
% setup the subdomains
%
DomainList.Ndomains = nxdomain*nydomain;
DomainList.dim = 2;
DomainList.LowerLeftCorner = zeros(DomainList.Ndomains,DomainList.dim);
DomainList.UpperRightCorner= zeros(DomainList.Ndomains,DomainList.dim);
DomainList.Bcs =  char(double('i')*ones(DomainList.Ndomains,DomainList.dim*2));
dx = 1/nxdomain;
dy = 1/nydomain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main idea here is to have smaller blocks near the center of the problem domain, as
% this is where the mesh is adaptively refined, and larger blocks near the edges where
% the mesh is coarse.
%fprintf('Using %d-times-%d subdomains with width reduction of 1/%d\n',nxdomain,nydomain,h);
%if rem(nxdomain,2) == 1 || rem(nydomain,2) == 1
%  error('number of domains must be even in x and y directions');
%end
%dx = zeros(nxdomain+1,1);
%halfdx = floor(nxdomain/2);
%dx(1) = -1;
%%dx(2) = -0.4;
%%for ii=3:halfdx
%%for ii=2:halfdx
%%  dx(ii) = dx(ii-1) + initialH / (h^(ii-2));
%%end
%for ii=halfdx:-1:2
%  dx(ii) = -1/h^(ii-1);
%end
%
%
%for ii=nxdomain+1:-1:halfdx+1
%  dx(ii) = -dx(nxdomain-ii+2);
%end
%
%dy = zeros(nydomain+1,1);
%halfdy = floor(nydomain/2);
%dy(1) = -1;
%%dy(2) = -0.4;
%%for ii=3:halfdy
%%for ii=2:halfdy
%%  dy(ii) = dy(ii-1) + initialH / (h^(ii-2));
%%end
%for ii=halfdy:-1:2
%  dy(ii) = -1/h^(ii-1);
%end
%
%for ii=nydomain+1:-1:halfdy+1
%  dy(ii) = -dy(nydomain-ii+2);
%end
%
%dx = (dx + 1)/2;
%dy = (dy + 1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i= 1:nxdomain
   for j= 1:nydomain
      DomainList.LowerLeftCorner((j-1)*nxdomain+i,1) = (i-1)*dx;
      DomainList.LowerLeftCorner((j-1)*nxdomain+i,2) = (j-1)*dy;
      DomainList.UpperRightCorner((j-1)*nxdomain+i,1)= (i  )*dx;
      DomainList.UpperRightCorner((j-1)*nxdomain+i,2)= (j  )*dy;
      %DomainList.LowerLeftCorner((j-1)*nxdomain+i,1) = dx(i);
      %DomainList.LowerLeftCorner((j-1)*nxdomain+i,2) = dy(j);
      %DomainList.UpperRightCorner((j-1)*nxdomain+i,1)= dx(i+1);
      %DomainList.UpperRightCorner((j-1)*nxdomain+i,2)= dy(j+1);
   end
end

%
% fix things up near the boundary. Make the upper right boundary a little
% bigger to make sure that nothing gets excluded when we do things like
% find(xx < UpperRightCorner(k,1). Also record the location of the boundarys
% in DomainList.Bcs. This is used to keep track of which domains are split
% among aggregates and which reside within an aggregate.
%
for j= 1:nydomain
      DomainList.Bcs((j-1)*nxdomain+1,1) = 'b';
      DomainList.UpperRightCorner(j*nxdomain,1)= 1.000001;
      DomainList.Bcs(j*nxdomain,2) = 'b';
end
for i= 1:nxdomain
   DomainList.Bcs(i,3) = 'b';
   DomainList.UpperRightCorner((nydomain-1)*nxdomain+i,2)= 1.000001;
   DomainList.Bcs((nydomain-1)*nxdomain+i,4) = 'b';
end

% Now set up the cross points by indicating the four subdomains that
% make up a cross point.
%%%% FIXME this looks bad....
CrossPointList.NCrossPoints =  (nxdomain-1)*(nydomain-1);
CrossPointList.Domains      =  zeros(CrossPointList.NCrossPoints,4);
for i= 1:nxdomain-1
   for j= 1:nydomain-1
      CrossPointList.Domains((j-1)*(nxdomain-1)+i,1)= (j-1)*nxdomain+i;
      CrossPointList.Domains((j-1)*(nxdomain-1)+i,2)= (j-1)*nxdomain+i+1;
      CrossPointList.Domains((j-1)*(nxdomain-1)+i,3)= (j  )*nxdomain+i;
      CrossPointList.Domains((j-1)*(nxdomain-1)+i,4)= (j  )*nxdomain+i+1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Read in data and massage a bit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(loadDataStr)
tic
eval(loadDataStr);
toc
fprintf('A: %d rows, %d cols, %d nnz\n',size(A,1),size(A,2),nnz(A));
fprintf('%d x %d domains for coarsening\n',nxdomain,nydomain);
fprintf('%d x %d domains for smoothing \n\n',smoonxdomain,smoonydomain);

if doom
  freeNode = setdiff(1:mesh.Nn,mesh.NBIDX);
end

xmax = max(mesh.node(:,1));
ymax = max(mesh.node(:,2));
xmin = min(mesh.node(:,1));
ymin = min(mesh.node(:,2));
xdelta = abs(xmax - xmin);
ydelta = abs(ymax - ymin);

%mesh.node = (mesh.node+1)/2;   % rescale so 0 <= coords <= 1
mesh.node(:,1) = (mesh.node(:,1)+abs(xmin))/xdelta;   % rescale so 0 <= xcoords <= 1
mesh.node(:,2) = (mesh.node(:,2)+abs(ymin))/ydelta;   % rescale so 0 <= ycoords <= 1

for i= 1:nxdomain+1
  newdx(i) = (i-1)*dx;
end
for i= 1:nydomain+1
  newdy(i) = (i-1)*dy;
end
dx = newdx; dy = newdy;

for i= 1:smoonxdomain+1
  newsmoodx(i) = (i-1)*smoodx;
end
for i= 1:smoonydomain+1
  newsmoody(i) = (i-1)*smoody;
end
smoodx = newsmoodx; smoody = newsmoody;

xx = mesh.node(:,1);
yy = mesh.node(:,2);
xx  = xx(freeNode);   % strip out boundary
yy  = yy(freeNode);

if size(A,1) ~= length(freeNode)
  %xx  = mesh.node(freeNode,1);   % strip out boundary
  %yy  = mesh.node(freeNode,2);
  A   = A(freeNode,freeNode);
end

if vizMesh
  figure
  gplot(A,[xx yy]);
  hold on;
  for ii=1:length(dx)
    oo = ones(length(dy),1);
    plot(dx(ii)*oo,dy,'r.-');
    hold on;
  end
  for ii=1:length(dy)
    oo = ones(length(dx),1);
    plot(dx,dy(ii)*oo,'r.-');
    hold on;
  end
  %%%%%%
  for ii=1:length(smoodx)
    oo = ones(length(smoody),1);
    plot(smoodx(ii)*oo,smoody,'c.--');
    hold on;
  end
  for ii=1:length(smoody)
    oo = ones(length(smoodx),1);
    plot(smoodx,smoody(ii)*oo,'c.--');
    hold on;
  end
end

if doom
  sol = rand(length(freeNode),1);
  rhs = A*sol;
  rhs = rhs / norm(rhs);
else
  rhs = A*mesh.solu(freeNode);
end

xx(xx == 0.) = 1e-7;     % put away from zero so some > 0 code in 
yy(yy == 0.) = 1e-7;     % the Agg and Pat factory get bc points

global special;
special = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Specify domain decomposition particulars
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mue_include;

AggFact= GeometricAggFact();
AggFact.SetGraphName('A');
ApproxPFact = TentativePFactoryEx([],AggFact,CoarseNSFactory());
EminSolver = CGEminSolver();
EminSolver.SetNumIterations(5);
Pfact  = EminPFactory(GeometricPatFact(), ConstraintFactory(), EminSolver, ApproxPFact, 'A');
PRfact = GenericPRFactory(Pfact,TransPFactory());
%pre    = TipSmoother('GaussSeidel', 2, 1.0);
%post   = TipSmoother('GaussSeidel', 2, 1.0);
pre    = Smoother('GaussSeidel', 2, 1.0);
post   = Smoother('GaussSeidel', 2, 1.0);
pre.SetBackwardSweep(false);
post.SetForwardSweep(false);

%
%  Construct and populate finest level with user information
%
Finest = Level(); %Finest.KeepAll(false);
Finest.Set('NullSpace', ones(size(A,1),1));
Finest.Set('A', Operator(A));
Finest.Set('xcoords',xx);
Finest.Set('ycoords',yy);
Finest.Set('DomainList',DomainList);
Finest.Set('CrossPointList',CrossPointList);
Finest.Set('Bcs',find( (xx<1.e-6) | (xx>1-1.e-6) |(yy<1.e-6) | (yy>1-1.e-6)));
%
% Set up a smoother corresponding to multiplicative Schwarz
%     Step 1: create nonoverlapping blocks
%     Step 2: grow blocks by 1 via blocks = spones( abs(A)*blocks);
%     Step 3: put block information into a subset used by Block Gauss-Seidel
%
n = length(xx);
Collection.NSubsets = smooDomainList.Ndomains;
onesA = spones(A);
tic
fprintf('Setting up smoother domains:   ');
for i=1:smooDomainList.Ndomains
   inside = ones(n,1);
   inside(xx <=  smooDomainList.LowerLeftCorner( i,1) ) = 0;
   inside(xx >   smooDomainList.UpperRightCorner(i,1) ) = 0;
   if smooDomainList.dim > 1,
      inside(yy <= smooDomainList.LowerLeftCorner( i,2) ) = 0;
      inside(yy >  smooDomainList.UpperRightCorner(i,2) ) = 0;
   end
   for jj=1:overlap
     inside = onesA*inside;
   end
   [aaa,bbb,ccc] = find(A(:,inside>0));
   Collection.Subsets(i)=CreateDOFSubset(-1,'Scattered',-1,-1,unique(aaa));
end
toc

% request for each TipLevelSmoother
Finest.Request('Tips');
Finest.Request('Tips');
Finest.Set('Tips',Collection);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);

status=MgHierarchy.FillHierarchy(PRfact, [], RAPFactory(),1,numDesiredLevels);

if numDesiredLevels > 1
  MgHierarchy.SetSmoothers(SmootherFactory(pre,post), 1, numDesiredLevels-1);
  A1 = MgHierarchy.Levels_{1}.Get('A');
  A2 = MgHierarchy.Levels_{2}.Get('A');
  fprintf('OC: (2-Level) = %4.2f\n', (nnz(A1.GetMatrixData)  + nnz(A2.GetMatrixData)) / nnz(A1.GetMatrixData) );
else
  MgHierarchy.SetCoarsestSolver(SmootherFactory(pre,post));
end

%MgHierarchy.TestMemoryLeak();
fprintf('CG/MG: maximum of %d iterations or stop tolerance of %g\n',maxIts,stopTol);
solTimer = tic;
[sol,flag,relres,iter,resvec]=pcg(A,rhs,stopTol,maxIts,...
   @(rhs)MgHierarchy.Iterate(rhs,1, zeros(n,1),ALLZEROS,1));
solTime = toc(solTimer);

fprintf('Nits = %d, ||r||/||rhs|| = %e,  %g seconds\n',length(resvec)-1,norm(rhs-A*sol)/norm(rhs),solTime);
if flag > 0, fprintf('CG error code %d: ',flag); end
switch flag

   case 0
     fprintf('Success! PCG converged.\n');
   case 1
     fprintf('maximum iterations hit without convergence\n');
   case 2
     fprintf('preconditioner is ill-conditioned\n');
   case 3
     fprintf('CG stagnated (two consecutive iterates are the same)\n');
   case 4
     fprintf('too large or small internal CG scalars\n');
end

%semilogy(resvec,linespec)
