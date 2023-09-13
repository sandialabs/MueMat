classdef UnitTest < handle
  % Little class to do unit testing.
  %
  % The main idea is everything gets passed in as strings and eval'd inside the tester.
  % This makes it easy to send in expressions to compare, auxiliary data, and comparision criteria.
  %

  properties (Access=private)
    numTests_
    label_
    method1_
    method2_
    dataString_
    evalString_
    userVerbLevel_
  end %properties

  methods
    function [this] = UnitTest()
      % Constructor.
      %
      % syntax: UnitTest();
      %
      this.numTests_ = 0;
      this.label_ = [];
      this.method1_ = [];
      this.method2_ = [];
      this.dataString_ = [];
      this.evalString_ = [];
      this.userVerbLevel_ = 0;
    end

    %%%%%%%

    function SetData(this,data)
      % Sets a string that when evaluated generates all auxiliary data necessary for tests.
      %
      % syntax: SetData(data)
      %
      % Example: ut.SetData('v = rand(50,1); v = v/norm(v);' ]);
      %
      this.dataString_ = data;
    end

    %%%%%%%

    function SetTest(this,label,method1,method2)
      % Registers a test to be run.
      %
      % syntax: SetTest(label,method1,method2)
      %
      % Example: ut.SetTest('Operator transpose', 'y1 = A''*v;', 'y2 = A.GetMatrixData()''*v;');
      %
      this.numTests_ = this.numTests_ + 1;
      this.label_{this.numTests_} = label;
      this.method1_{this.numTests_} = method1;
      this.method2_{this.numTests_} = method2;
    end

    %%%%%%%

    function [method1,method2] = GetTest(this,ii)
      if ii < 0 || ii > this.numTests_
        error('No such test');
      end
      method1 = this.method1_{ii};
      method2 = this.method2_{ii};
    end

    %%%%%%%

    function SetOutputLevel(this,outputLevel)
      % Sets amount of output the user wants (0,1, or 2)
      %
      % syntax:  SetOutputLevel(outputLevel)
      %
      this.userVerbLevel_ = outputLevel;
    end

    %%%%%%%

    function SetEvaluator(this,evalString)
      % Sets a string that when evaluated determines whether the test passed.
      % All tests will use use this criterion.
      %
      % syntax: SetEvaluator(evalString)
      %
      % Example:  ut.SetEvaluator('norm(r) < 1e-12');
      %
      this.evalString_ = evalString;
    end

    %%%%%%%

    function Run(this)
      % Executes all registered tests.
      %
      % syntax: Run()
      %
      fprintf('Running %d tests\n',this.numTests_);
      eval(this.dataString_);
      for ii=1:this.numTests_
        fprintf('%3d /%3d %-40s',ii,this.numTests_,this.label_{ii});
        eval(this.method1_{ii});
        eval(this.method2_{ii});
        %if (norm(y1-y2) < 1e-10) TEST_PASSED=true;
        if ( eval(this.evalString_) ) TEST_PASSED=true;
        else                     TEST_PASSED=false; end
        VERBOSE = this.DetermineOutputLevel(TEST_PASSED,this.userVerbLevel_);
        if VERBOSE
          fprintf('\n\n%s\n',this.method1_{ii});
          disp(y1)
          fprintf('\n\n%s\n',this.method2_{ii});
          disp(y2)
        end
        if TEST_PASSED fprintf(' ... passed\n');
        else           fprintf(' ... FAILED (%g)\n',norm(y1-y2)); end
        if VERBOSE && ~TEST_PASSED, keyboard; end
      end %for ii
    end %Run()

  end %public methods

  methods (Access=private)

    %%%%%%%

    function [finalVerbLevel] = DetermineOutputLevel(this,testPassed,userVerbLevel)
      finalVerbLevel = (~testPassed && userVerbLevel>0) || (userVerbLevel>1);
    end % CheckVerbose()

  end %private methods

end %classdef UnitTest
