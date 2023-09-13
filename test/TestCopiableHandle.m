% Test of CopiableHandle Copy() method.

function TestCopiableHandle()
    SetHomeDir;
    P=path; path([MUEMAT_ROOT_DIR,'/test/TestCopiableHandle'],P);

    try
      Tests();
    catch Exception
      fprintf('%s\n',Exception.message);
    end

    path(P);
end

function Tests()

%% Basic tests

% Create an instance of Foo which stores 4 integers
F1 = Foo(1,2);
F1.SetBaseA(3);
F1.SetBaseB(4);

% Create copies of F1
F2 = Foo(F1);
F3 = F1.Copy();
F4 = F1;

% Check handle vs deep copies
unitTest('F1 ~= F2');
unitTest('F1 ~= F3');
unitTest('F1 == F4');

% Check if all properties values have been copied
checkFooValue(F1, F2, '==');
checkFooValue(F1, F3, '==');

clear;

%% Recursive tests

% Create an instance of Foo which stores 4 Foo objects
F1 = Foo();
F1.SetFooA(Foo(1,1));
F1.SetFooB(Foo(2,2));
F1.SetBaseA(Foo(3,3));
F1.SetBaseB(Foo(4,4));

% Create copies of F1
F2 = Foo(F1);
F3 = F1.Copy();
F4 = F1;

% Check handle vs deep copies of F1 properties (== F1.FooA_/F1.FooB_/F1.BaseA_/F1.BaseB_)
checkFooValue(F1, F2, '~=');
checkFooValue(F1, F3, '~=');
checkFooValue(F1, F4, '==');

% Check if all properties values of objects F1.FooA_/F1.FooB_/F1.BaseA_/F1.BaseB_ have been copied
checkFooValue(F1.GetFooA(),  F2.GetFooA(),  '==');
checkFooValue(F1.GetFooB(),  F2.GetFooB(),  '==');
checkFooValue(F1.GetBaseA(), F2.GetBaseA(), '==');
checkFooValue(F1.GetBaseB(), F2.GetBaseB(), '==');
%
checkFooValue(F1.GetFooA(),  F3.GetFooA(),  '==');
checkFooValue(F1.GetFooB(),  F3.GetFooB(),  '==');
checkFooValue(F1.GetBaseA(), F3.GetBaseA(), '==');
checkFooValue(F1.GetBaseB(), F3.GetBaseB(), '==');

fprintf('Test passed !\n');

end

function checkFooValue(F1, F2, op)
  unitTest(['F1.GetFooA()  ' op ' F2.GetFooA()']);
  unitTest(['F1.GetFooB()  ' op ' F2.GetFooB()']);
  unitTest(['F1.GetBaseA() ' op ' F2.GetBaseA()']);
  unitTest(['F1.GetBaseB() ' op ' F2.GetBaseB()']);
end

function unitTest(string)
  bool = evalin('caller',string);
  if bool ~= 1
     error(['unitTest error: ' string]);
  end
end  

% TODO: Add test for SetBaseA(vector + cell)