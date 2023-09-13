function OperatorTest(nrows,VERBOSE_TEST)
  %%
  % syntax: OperatorTest(nrows,VERBOSE_TEST)
  %
  %   nrows        -- number of rows in test matrix       (default=10)
  %   VERBOSE_TEST -- controls amount of screen output    (default=0)
  %                   0 => minimal output, only report failing tests
  %                   1 => minimal output, stop in failing tests
  %                   2 => details about all tests
  %
  % Unit tests to check that various pieces of the Operator class work properly.
  % Testing is based on performing the operation in two different ways on the same
  % random vector and comparing the norms of the two results.
  %
  % See also: UnitTest

  if exist('nrows') ~= 1 || isempty(nrows), nrows = 10;       end
  if exist('VERBOSE_TEST') ~= 1,            VERBOSE_TEST = 0; end

  density = 0.75;
  % data used in the tests
  %A = BuildRandomOp(nrows,-10,10,density);
  %v = rand(nrows,1);
  %v = v/norm(v);
  %B = BuildRandomOp(nrows,-nrows,nrows,density);
  %alpha = round(rand*10)+1;
  dataString = ['A = BuildRandomOp(' num2str(nrows) ',-10,10,' num2str(density) ');' ];
  dataString = [ dataString 'v = rand(' num2str(nrows) ',1);' ];
  dataString = [ dataString 'v = v/norm(v);' ];
  dataString = [ dataString 'B = BuildRandomOp(' num2str(nrows) ',-10,10,' num2str(density) ');' ];
  dataString = [ dataString 'alpha = round(rand*10)+1;' ];

  evalString = 'norm(y1-y2) < 1e-10';

  ut = UnitTest();
  ut.SetData(dataString);
  ut.SetEvaluator(evalString);
  ut.SetOutputLevel(VERBOSE_TEST);

  ut.SetTest('Operator vector multiplication', ...
             'y1 = A*v;', ...
             'y2 = A.GetMatrixData() * v;');

  ut.SetTest('Operator/Operator addition', 'C=A+B; y1 = C*v;', 'y2 = A*v + B*v;');

  ut.SetTest('Operator/Matrix addition', 'C=A+B.GetMatrixData(); y1 = C*v;', 'y2 = A*v + B*v;');

  ut.SetTest('Operator/Operator subtraction', 'C=A-B; y1 = C*v;', 'y2 = A*v - B*v;');

  ut.SetTest('Operator scaling', 'C=A.Copy(); C.Scale(alpha); y1 = C*v;', 'y2 = alpha*(A*v);');

  ut.SetTest('Operator-operator multiplication', 'y1 = (A*B)*v;', 'y2 = A*(B*v);');

  ut.SetTest('Operator transpose', 'y1 = A''*v;', 'y2 = A.GetMatrixData()''*v;');

  ut.SetTest('Extract Operator diagonal', 'tt=GetDiagonal(A); y1 = tt.GetMatrixData();', 'y2 = diag(A.GetMatrixData());');

  ut.SetTest('Apply Operator diagonal', 'tt=GetDiagonal(A); FactorBlkDiag(tt); y1 = tt \ v;', 'y2 = diag(diag(A.GetMatrixData())) \ v;');

  ut.SetTest('Left divide (v is a vector)', 'y1 = A \ v;', 'y2 = A.GetMatrixData() \ v;');

  ut.SetTest('Left divide (v is a MultiVector)', 'mv = MultiVector([],[],v); y1 = A \ mv;', 'y2 = A.GetMatrixData() \ v;');

  ut.SetTest('Left divide (v is an Operator)', 'mat = Operator(v); y1 = A \ mat; y1 = y1.GetMatrixData();', 'y2 = A.GetMatrixData() \ v;');

  ut.SetTest('Matrix MultiVector multiplication', 'y1 = A*v;', 'y2 = A.GetMatrixData() * v;');

  ut.Run();

end % OperatorTest()
