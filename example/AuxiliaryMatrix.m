%
% Test driver to demonstrate projection capability for auxiliary matrices and vectors.
% Two different Laplace problems are used (constant block size and
% variable block size). Currently, MueMat is only able to project
% coordinates, and an auxiliary matrix. Read the comments associated with
% Level.ProjectInterface() to see how other matrices/vectors
% can be added to the list of projected items.
%

srand;

% Set things that are held constant through all tests
AmalgamateDropFact = CoalesceDropFactory();
Pfact              = SaPFactory(TentativePFactory(AmalgamateDropFact,AggregationFactory()));

Rfact        = TransPFactory();
PRfact       = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(3);
Acfact       = RAPexFactory();
NBlks        = 81;
coords       = (1:NBlks)'/NBlks;

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Variable Block Size Problem projecting coordinate values and\n');
fprintf('an auxiliary matrix which happens to be equal to Amat\n');

Bsizes =  -sort(-ceil(6*rand(NBlks,1))); Bsizes(Bsizes < 3) = 3;
VarBlkPtr(1) = 1;
for i=2:NBlks+1, VarBlkPtr(i) = VarBlkPtr(i-1) + Bsizes(i-1); end

Amat   = BuildLaplace1DBlk(-1,VarBlkPtr,NBlks);
Finest = Level();
Finest.KeepAll(false);

Finest.Set('A', Amat);

Finest.Keep('NullSpace');
Finest.Set('NullSpace', BuildNullSpace(Amat));

Finest.Keep('xcoords');
Finest.Set('xcoords', coords);

Finest.Keep('AuxMatrix');
Finest.Set('AuxMatrix', Amat);

AmalgamateDropFact.SetAName('AuxMatrix');

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest.Copy(),1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 5);

fprintf('Hierarchy filled. Printing some projected coordinates\n');
full(MgHierarchy.Levels_{3}.Get('xcoords'))',
fprintf('\nNow comparing project auxiliary matrix with coarse Amat\n');
fprintf('nnz( Amatcoarse - projected Aux ) = %d\n',...
nnz(MgHierarchy.Levels_{4}.Get('A').GetMatrixData()-...
MgHierarchy.Levels_{4}.Get('AuxMatrix').GetMatrixData()));
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Variable Block Size Problem projecting coordinate values and\n');
fprintf('an auxiliary matrix which has 1 dof/node (built by coalescing.\n');
fprintf('entries in Amat).\n');
ExtractFact = CoalesceDropFactory();
ExtractFact.SetPostDropSpecifications([],[]);
ExtractFact.SetPreDropSpecifications([],[]);
ExtractFact.SetAmalgSpecifications(@PickFirst,[]);
ExtractFact.SetNonBinary();

Finest.Delete('AuxMatrix'); % avoid a warning
Finest.Keep('AuxMatrix');
Finest.Set('AuxMatrix', ExtractFact.Build_(Amat));
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 5);
fprintf('\nPrinting a projected auxiliary matrix. This would look like a\n');
fprintf('scaled Laplace problem if the tenative prolongator was used for\n');
fprintf('projecting auxiliary information. Unfortunately, there is no easy\n');
fprintf('way to currently turn-off prolongator smoothing. Another problem\n');
fprintf('is that the RAP factory works with Generic Data Buckets ... so\n');
fprintf('we need to do something smart if it is to somehow use a tentative\n');
fprintf('prolongator. There might be a way to have it query the databucket\n');
fprintf('and get the grid transfer operators out?\n');

full(MgHierarchy.Levels_{3}.Get('AuxMatrix').GetMatrixData()),

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Constant Block Size Problem projecting coordinate values and\n');
fprintf('an auxiliary matrix which happens to be equal to Amat\n');
ConstBlkSize= 4; VarBlkPtr = [];
Amat   = BuildLaplace1DBlk(ConstBlkSize,VarBlkPtr,NBlks);
Finest = Level();
Finest.KeepAll(false);
Finest.Set('A', Amat);
Finest.Keep('NullSpace');
Finest.Set('NullSpace', BuildNullSpace(Amat));
Finest.Keep('xcoords');
Finest.Set('xcoords', coords);
Finest.Keep('AuxMatrix');
Finest.Set('AuxMatrix', Amat);
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest.Copy(),1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 5);

fprintf('Hierarchy filled. Printing some projected coordinates\n');
full(MgHierarchy.Levels_{3}.Get('xcoords'))',
fprintf('\nNow comparing project auxiliary matrix with coarse Amat\n');
fprintf('nnz( Amatcoarse - projected Aux ) = %d\n',...
nnz(MgHierarchy.Levels_{4}.Get('A').GetMatrixData()-...
MgHierarchy.Levels_{4}.Get('AuxMatrix').GetMatrixData()));
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Constant Block Size Problem projecting coordinate values and\n');
fprintf('an auxiliary matrix which has 1 dof/node (built by coalescing.\n');
fprintf('entries in Amat).\n');
ExtractFact = CoalesceDropFactory();
ExtractFact.SetPostDropSpecifications([],[]);
ExtractFact.SetPreDropSpecifications([],[]);
ExtractFact.SetAmalgSpecifications(@PickFirst,[]);
ExtractFact.SetNonBinary();
Finest.Delete('AuxMatrix'); % avoid a warning
Finest.Keep('AuxMatrix');
Finest.Set('AuxMatrix', ExtractFact.Build_(Amat));
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 5);
fprintf('\nPrinting a projected auxiliary matrix. This would look like a\n');
fprintf('scaled Laplace problem if the tenative prolongator was used for\n');
fprintf('projecting auxiliary information. Unfortunately, there is no easy\n');
fprintf('way to currently turn-off prolongator smoothing. Another problem\n');
fprintf('is that the RAP factory works with Generic Data Buckets ... so\n');
fprintf('we need to do something smart if it is to somehow use a tentative\n');
fprintf('prolongator. There might be a way to have it query the databucket\n');
fprintf('and get the grid transfer operators out?\n');
full(MgHierarchy.Levels_{3}.Get('AuxMatrix').GetMatrixData()),
