%%% The operator and multivector class tutorial
% This tutorial presents the basic usage of Operator and Vector objects.

%%
% See also: <ref_Operator.html Operator>, <ref_MultiVector.html MultiVector>

clear;

%% Build a scalar operator
Amat = gallery('poisson',50); % a matlab matrix
opScalar = Operator(Amat);

disp(Amat - opScalar.GetMatrixData());

%% Build a 2x2-block operator
op2x2Blk = Operator(Amat,2,2);

if op2x2Blk.GetRowMap().HasConstBlkSize()
  fprintf('Stats of op2x2Blk:\n');
  fprintf('  nDOFs      : %d\n', op2x2Blk.GetRowMap().NDOFs());
  fprintf('  nBlk       : %d\n', op2x2Blk.GetRowMap().NNodes()); % number of block per row
  fprintf('  blkSize    : %d\n', op2x2Blk.GetRowMap().ConstBlkSize());
end

%% Build a 2x2-block operator using a pre-defined Map object
ConstBlkSize = 2;                            % block size
nBlk         = size(Amat,1) / ConstBlkSize;  % number of blocks
ConstBlkMap  = Map(nBlk, ConstBlkSize);
opConstBlk   = Operator(Amat, ConstBlkMap, ConstBlkMap, @MatlabApply);

%% Build a variable block size operator
% On this example, the size of each blocks is defined randomely.
%
% See also: <matlab:doc('BuildLaplace1DBlk') BuildLaplace1DBlk>
srand;                % 'srand' initializes the random number generator to ensure reproducible results
r = size(Amat,1);
blkSize = zeros(1,r); % array preallocation
i = 1;                % i: current block number
while r ~= 0
  blkSize(i) = randi([1,min(10,r)],1); % Generate random block sizes, without exceeding the matrix size.
  r = r - blkSize(i);
  i = i + 1;
end
nBlk = i;
blkSize = blkSize(1:nBlk); % truncate according to the number of blocks.

VarBlkMap  = Map(nBlk, blkSize);
opVarBlk   = Operator(Amat, VarBlkMap, VarBlkMap, @MatlabApply);

if opVarBlk.GetRowMap().HasVariableBlkSize()
  fprintf('Stats of opVarBlk:\n');
  fprintf('  nDOFs      : %d\n', opVarBlk.GetRowMap().NDOFs());
  fprintf('  nBlk       : %d\n', opVarBlk.GetRowMap().NNodes());
  fprintf('  maxBlkSize : %d\n', opVarBlk.GetRowMap().MaxBlkSize());
end

%% Build vectors

% 1 dense vector of length 50 and initialized to zero
u = MultiVector(50,1);

% 1 dense vector of length 50 and initialized from user data
m = zeros(1,50); % a matlab array
v = MultiVector([],[],m);

% 5 dense vector of length 50 and initialized from user data
m = zeros(5,50); % a matlab array
w = MultiVector([],[],m);

fprintf('w.NumVectors() : %d\n',w.NumVectors());

%% Copy vectors or operators
u = MultiVector(50,1); % u = [0 0 0 ... 0]

% Handle copy
v = u;  % r and w refers to the same data in memory
fprintf('norm(v) : %g \t norm(u) : %g\n', u.norm(1), v.norm(1));

v.PutScalar(1); % v = [1 1 1 ... 1]
fprintf('norm(v) : %g \t norm(u) : %g\n', u.norm(1), v.norm(1));

% Data copy
v = MultiVector(u); % v is a copy of u
w = u.Copy();       % w is a copy of u (using another syntax)

v.PutScalar(2); % v = [1 1 1 ... 1]
fprintf('norm(v) : %g \t norm(u) : %g\n', u.norm(1), v.norm(1));

%% Operator views
% User can access the matrix as if it were stored in different
% formats, independent of the actual underlying storage format.
% Eachway of accessing the matrix is called a "view".
% A view allows for associating with the matrix a particular row map,
% column map, matrix-vector multiply, and notion of matrix diagonal.

% Build a 5x5-block operator
op = Operator(Amat,5,5);

fprintf('Stats of op:\n');
fprintf('  current view    :     ''%s'' (blkSize=%d)\n', op.CurrentView(), op.GetRowMap().ConstBlkSize());
fprintf('  available views : '); disp(op.GetViewList())

% Add a view, using a Map defined previously
op.CreateView('2x2',ConstBlkMap, ConstBlkMap, @MatlabApply);

fprintf('  available views : '); disp(op.GetViewList())

% Switch to the new view '2x2'
% See also: SwitchToDefaultView() and SwitchToPointView();
oldview = op.SwitchToView('2x2'); % output 'oldview' is optional

fprintf('  current view    :     ''%s'' (blkSize=%d)\n', op.CurrentView(), op.GetRowMap().ConstBlkSize());
fprintf('  old view        :     ''%s''\n', oldview);

%% Operator diagonals
% Get the scalar diagonal
op.SwitchToPointView();
opDiag = op.GetDiagonal(); % opDiag is a single view operator

disp(diag(Amat) - opDiag.GetMatrixData());

% Get the block diagonal
op.SwitchToDefaultView();
op.SwitchToView('2x2');
opBlkDiag = op.GetDiagonal(); % opBlkDiag is a single view operator

disp(size(opBlkDiag.GetMatrixData())); % diagonal block (2x2) are stored continuously

%% Loops 'for each operator blocks'
% This example creates a cell array where each cell correspond to a
% block of the inital matrix.
%
op   = Operator(Amat,5,3);  % note that RowMap and ColMap can be different.
Amat = op.GetMatrixData();

%
nRowBlk = op.GetRowMap().NNodes();
nColBlk = op.GetRowMap().NNodes();
Blk  = cell(nRowBlk,nColBlk);

if op.GetRowMap().HasConstBlkSize()

  RowBlkSize = op.GetRowMap().ConstBlkSize();
  ColBlkSize = op.GetColMap().ConstBlkSize();

  for i=1:nRowBlk
    fRow = (i-1)*RowBlkSize+1;
    lRow = fRow + RowBlkSize-1;

    for j=1:nColBlk
      fCol = (j-1)*ColBlkSize+1;
      lCol = fCol + ColBlkSize-1;

      c{i,j} = Amat(fRow:lRow, fCol:lCol);
    end % for j
  end % for i

end % if

% Another example with a variable block operator
op  = Operator(Amat, VarBlkMap, VarBlkMap, @MatlabApply);

nRowBlk = op.GetRowMap().NNodes();
nColBlk = op.GetRowMap().NNodes();
Blk = cell(nRowBlk,nColBlk);

if op.GetRowMap().HasVariableBlkSize()

  VarRowBlkPtr = op.GetRowMap().Vptr();
  VarColBlkPtr = op.GetColMap().Vptr();
  for i=1:nRowBlk
    fRow = VarRowBlkPtr(i);
    lRow = VarRowBlkPtr(i+1)-1;

    for j=1:nColBlk
      fCol = VarColBlkPtr(j);
      lCol = VarColBlkPtr(j+1)-1;

      c{i,j} = Amat(fRow:lRow, fCol:lCol);
    end % for j
  end % for i

end % if