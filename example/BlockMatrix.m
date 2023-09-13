%
%
%
% Something to exercise the AMG code. I use random numbers to try out
% different things. In particular,
%
%     z = rand(4,1);
%
%     z(1)*100  Gives the block dimension of the system to be solved.
%
%     z(2)      If >= .5, a constant blocksized problem is built.
%               Othersie, a variable blocksized problem is constructed.
%               The size of individual blocks range between 9 and 3
%               which are sorted from highests to lowest sizes
%               (BuildLaplace1DBlk.m used to require this. I can't remember
%               if it still does).
%               Note: one has to be a little careful with having
%                     aggregates big enough so that the QR works.
%                     If any 9x9 blocks are present in the system,
%                     the dimension of the null space will be
%                     9. Since some blocks can be as small as 3x3,
%                     aggregates should have at least 3 nodes
%                     so that the smallest matrix given to QR is
%                     a 9x9 matrix (9 = 3x3).
%
%     z(3)*10   Gives the block size of constant blockedsized
%               problems (i.e. when z(2) >= .5).
%
%     z(4)      If < .5, arbitrary/random blocks are constructed
%               for the block diagonal and used with block
%               Gauss-Seidel.
%               Otherwise, block Gauss-Seidel uses the standard
%               block diagonal matrix of the discretization
%               operator.


for kkk=1:5   % run 10 problems.
  clear
  mue_include
  srand;
  SolStatus = NOTALLZEROS;

  % Build the problem matrix
  VarBlkPtr = [];

  z=rand(4,1); %z = [1;.2;1; 1];
  NBlks       = ceil(200*z(1)) + 3;

  if z(2) < .5,
     ConstBlkSize=-1;  Bsizes =  -sort(-ceil(9*rand(NBlks,1)));
     Bsizes(Bsizes < 3) = 3;
     VarBlkPtr(1) = 1;
     for i=2:NBlks+1, VarBlkPtr(i) = VarBlkPtr(i-1) + Bsizes(i-1); end
  else
     ConstBlkSize= ceil(10*z(3));
  end
  Amat   = BuildLaplace1DBlk(ConstBlkSize,VarBlkPtr,NBlks);

  % Set options
  numDesiredLevels = 4;         % number of AMG levels

  Pfact             = SaPFactory(TentativePFactory());
  Rfact             = TransPFactory();
  Acfact            = RAPFactory();
  Smoo              = Smoother('GaussSeidel', 2, 1);
  if z(4) < .5, Smoo.SetDiagonalView('Random NonOverlapping'); end
  GSFactory         = SmootherFactory(Smoo);

  %
  %  Construct and populate finest level with user information
  %
  Finest = Level();
  Finest.KeepAll(false);
  Finest.Set('A', Amat);
  Finest.Set('NullSpace', BuildNullSpace(Amat));


  MgHierarchy = Hierarchy();
  MgHierarchy.SetOutputLevel(1);
  MgHierarchy.SetLevel(Finest,1);
  MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
  MgHierarchy.SetSmoothers(GSFactory);

  rhs = rand(Amat.GetRowMap().NDOFs(),1);
  sol = zeros(Amat.GetRowMap().NDOFs(),1);  SolStatus = ALLZEROS;
  sol = MgHierarchy.Iterate(rhs, 9, sol, SolStatus);
end
