% Test setup phase reuse capability of SmootherFactory, Hybrid2x2Smoother and MergedSmoother

function SmoothingTestReUse()
  global count; count=0;
  
  global Level; Level=[];
  
  fprintf('SmoothingTestReUse: \n');
  
  TestSmootherFactory();
  TestTestSmoother();
  TestHybrid2x2Smoother(@BuildHybridSmootherFactory);
  TestHybrid2x2Smoother(@BuildOldHybridSmootherFactory);
  TestHybrid2x2Smoother(@BuildMergedSmootherFactory);
  
  fprintf('Number of passed tests: %d\n', count);
end

function TestSmootherFactory()
  TypeA  = TestSmoother('TypeA');
  TypeA_ = TestSmoother('TypeA'); % another one
  TypeB  = TestSmoother('TypeB');
  
  % TEST : PreSmoo != PostSmoo (-> no reuse)
  runTest(TypeA, TypeB, 2);
  
  % TEST : PreSmoo == PostSmoo (-> reuse)
  runTest(TypeA, TypeA, 1);  % test with two pointers to same smoother prototype
  runTest(TypeA, TypeA_, 1); % test with two distinct but equal smoother prototype
end

function TestTestSmoother()
  TypeA1  = TestSmoother('TypeA', 1);
  TypeA2  = TestSmoother('TypeA', 2);
  
  % TEST : PreSmoo == PostSmoo but != parameters (-> no reuse)
  runTest(TypeA1, TypeA2, 2);
end

function TestHybrid2x2Smoother(func)
  TypeA  = TestSmoother('TypeA');
  TypeA_ = TestSmoother('TypeA'); % another one
  TypeB  = TestSmoother('TypeB');
  TypeB_ = TestSmoother('TypeB'); % another one
  TypeC  = TestSmoother('TypeC');
  TypeD  = TestSmoother('TypeD');
  
  % TEST: PreSmoo != PostSmoo (-> no reuse)
  runTestHybrid(func,TypeA, TypeB, TypeC, TypeD, 4);
  
  % TEST: PreSmoo == PostSmoo (-> reuse)
  runTestHybrid(func,TypeA, TypeB, TypeA, TypeB, 2);
  runTestHybrid(func,TypeA, TypeB, TypeA_, TypeB_, 2);
  
  % TEST: (-> partial reuse)
  runTestHybrid(func,TypeA, TypeB, TypeA, TypeD, 3);  % PreSmooOne == PostSmooOne
  runTestHybrid(func,TypeA, TypeB, TypeA_, TypeD, 3); % idem
  runTestHybrid(func,TypeA, TypeB, TypeC, TypeB, 3);  % PreSmooTwo == PostSmooTwo
  runTestHybrid(func,TypeA, TypeB, TypeC, TypeB_, 3); % idem
  
  % TEST PreOne == PostTwo and PostOne == PreTwo -> no reuse
  runTestHybrid(func,TypeA, TypeB, TypeC, TypeA, 4); % only PreOne == PostTwo
  runTestHybrid(func,TypeA, TypeB, TypeB, TypeA, 4);
  runTestHybrid(func,TypeA, TypeB, TypeB, TypeC, 4); % only PostOne == PreTwo
  
  % with != parameters:
  TypeA1 =  TestSmoother('TypeA', 1);
  TypeA2 =  TestSmoother('TypeA', 2);
  TypeB1 =  TestSmoother('TypeB', 1);
  TypeB2 =  TestSmoother('TypeB', 2);
  
  % TEST: PreSmoo == PostSmoo but != parameters (-> no reuse)
  runTestHybrid(func,TypeA1, TypeB1, TypeA2, TypeB2, 4);
  
  % TEST: (-> no partial reuse)
  runTestHybrid(func,TypeA1, TypeB, TypeA2, TypeD,  4); % PreSmooOne == PostSmooOne but != parameters
  runTestHybrid(func,TypeA,  TypeB, TypeC,  TypeB2, 4); % PreSmooTwo == PostSmooTwo but != parameters
end

function runTest(PreSmooProto, PostSmooProto, NumOfSetupPhase)
  global count; count=count+1;
  
  SmooFact  = SmootherFactory(PreSmooProto, PostSmooProto);
  [PreSmoo, PostSmoo] = SmooFact.Build([],[]);
  
  if (PreSmooProto.SetupDone() + PostSmooProto.SetupDone() ~= 0) , error('Test'); end
  
  if (PreSmoo.SetupDone() + PostSmoo.SetupDone() ~= NumOfSetupPhase), error('Test'); end
end

function runTestHybrid(func, PreSmooProtoOne, PreSmooProtoTwo, PostSmooProtoOne, PostSmooProtoTwo, NumOfSetupPhase)
  global count; count=count+1;
  
  SmooFact = func(PreSmooProtoOne, PreSmooProtoTwo, PostSmooProtoOne, PostSmooProtoTwo);
  Level    = BuildLevel();
  
  [PreSmoo, PostSmoo] = SmooFact.Build(Level,[]);
  
  if (PreSmoo.GetSmootherOne().isSetup() + PreSmoo.GetSmootherTwo.isSetup() + PostSmoo.GetSmootherOne().isSetup() + PostSmoo.GetSmootherTwo.isSetup() ~= 4), error('Test'); end

  %[PreSmoo.GetSmootherOne().SetupDone(), PreSmoo.GetSmootherTwo.SetupDone(), PostSmoo.GetSmootherOne().SetupDone(), PostSmoo.GetSmootherTwo.SetupDone()]
  if (PreSmoo.GetSmootherOne().SetupDone() + PreSmoo.GetSmootherTwo.SetupDone() + PostSmoo.GetSmootherOne().SetupDone() + PostSmoo.GetSmootherTwo.SetupDone() ~= NumOfSetupPhase), error('Test'); end
end

function [lvl] = BuildLevel()
  global Level;
  
  if isempty(Level)
    
    n      = 30;
    Amat   = BuildLaplace2D(n);
    Level = Level();
    Level.Set('A', Amat);
    Level.Set('NullSpace', BuildNullSpace(Amat));
    
    NumR = 50;
    Bsize = 1;
    AmatData=Amat.GetMatrixData();
    AOne = Operator(AmatData(1:NumR,1:NumR),Bsize,Bsize);
    Level.Set('Arr', AOne);
    
  end
  
  lvl = Level;
end

function [SmooFact] = BuildHybridSmootherFactory(PreSmooProtoOne, PreSmooProtoTwo, PostSmooProtoOne, PostSmooProtoTwo)
  nMainIts = 5;
  MiddleOne = 5;
  StartTwo = 5;
  MiddleTwo = 5;
  EndTwo = 5;
  
  PreSmooProto = Hybrid2x2Smoother(nMainIts, MiddleOne, ...
    StartTwo, MiddleTwo, EndTwo, ...
    PreSmooProtoOne, PreSmooProtoTwo);
  
  PostSmooProto = Hybrid2x2Smoother(nMainIts, MiddleOne, ...
    StartTwo, MiddleTwo, EndTwo, ...
    PostSmooProtoOne, PostSmooProtoTwo);
  
  SmooFact  = SmootherFactory(PreSmooProto, PostSmooProto);
end

function [SmooFact] = BuildOldHybridSmootherFactory(PreSmooProtoOne, PreSmooProtoTwo, PostSmooProtoOne, PostSmooProtoTwo)
  nMainIts = 5;
  MiddleOne = 5;
  StartTwo = 5;
  MiddleTwo = 5;
  EndTwo = 5;
  
  PreSmooProto = OldHybrid2x2Smoother(nMainIts, MiddleOne, ...
    StartTwo, MiddleTwo, EndTwo, ...
    PreSmooProtoOne, PreSmooProtoTwo);
  
  PostSmooProto = OldHybrid2x2Smoother(nMainIts, MiddleOne, ...
    StartTwo, MiddleTwo, EndTwo, ...
    PostSmooProtoOne, PostSmooProtoTwo);
  
  SmooFact  = OldHybrid2x2SmootherFactory(PreSmooProto, PostSmooProto);
end

function [SmooFact] = BuildMergedSmootherFactory(PreSmooProtoOne, PreSmooProtoTwo, PostSmooProtoOne, PostSmooProtoTwo)
  PreSmooProto = MergedSmoother(PreSmooProtoOne, PreSmooProtoTwo);
  PostSmooProto = MergedSmoother(PostSmooProtoOne, PostSmooProtoTwo);
  SmooFact  = SmootherFactory(PreSmooProto, PostSmooProto);
end


