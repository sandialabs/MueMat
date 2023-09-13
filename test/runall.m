%
% Function that gives a convenient way to run lots of tests in sequence.
%
function runall(varargin)
  tests = {};
  
  tests = {   'TestPGAMG' ...
              'SIMPLESmooTest' ...
              'ConvDiff1D' ...
              'DropTest' ...
              'GeoDDtest' ...
              'NonSym' ...
              'MergedSmootherWithSpecialBlkGS' ...
              'ImplicitSchur' ...
              'ExplicitSchur' ...
              'CrossFactorySpec' ...
              'SmoothedAggregation' ...
              'BlockMatrix' ...
              'Laplace' ...
              'Simple' ...
              'TestSmoothing' ...
              'TestSmoothingILU' ...
              'TestSmoothingChebychev' ...
              'Views' ...
              'AuxiliaryMatrix' ...
              'ElasticityTest(2)' ...
              'ElasticityTest(3)' ...
              'TransferOpTest' ...
              'MultiVectorTest' ...
              'TestReUseOptions'
             };

  if nargin > 0, tests = addTests(tests, varargin{:}); end
  tests = unique(tests);
  numTests = length(tests);

  errors = {};

  for ii=1:length(tests)
    fprintf('\n=================================================================================\n');
    fprintf('Running example %2d of %2d:  %s\n',ii,numTests,tests{ii});
    fprintf('=================================================================================\n');

    try
      invokeExample(tests{ii});
    catch Exception
      fprintf('\nExample %2d:  %s: ERROR !\n',ii,tests{ii});
      fprintf('Exception:%s\n',Exception.identifier);
      errors{length(errors)+1} = tests{ii};
    end
  end
  
  fprintf('\n\n\n');
  fprintf('=================================================================================\n');
  fprintf('Summary:\n');
  fprintf('=================================================================================\n');  
      
  if (isempty(errors)) 
    fprintf('All tests passed.\n');
  else 
    fprintf('Failed %d test(s) of %d:\n', length(errors), length(tests));
    for ii=1:length(errors)
      fprintf('\t- %s\n',errors{ii});
    end
  end
  
  fprintf('=================================================================================\n');
  fprintf('=================================================================================\n');
end

function invokeExample(example)
  eval(example);
end

function tests = addTests(tests, varargin)
  SetHomeDir
  
  for d=1:size(varargin,1)
    filelist = dir([MUEMAT_ROOT_DIR '/' varargin{d} '/*.m']);
    for i=1:size(filelist,1)
      if ~filelist(i).isdir && ~strcmp(filelist(i).name, 'runall.m')
        tests{length(tests)+1} = filelist(i).name(1:end-2);
      end
    end  
  end

end
