%TODO: review what 'numDesiredLevel' means
% loops "for i=firstLevel::numDesiredLevel" in FillHierarchy() suggest that
% numDesiredLevel is mix up with the last level ID.
% If firstLevel == 2 and numDesiredLevel == 3, then lastLevel should be
% == 4 (2+3-1)

classdef Hierarchy < VerboseObject
  %Hierarchy Provides methods to build an multigrid hierarchy and apply multigrid cycles
  %
  % Allows users to manually populate operators at different levels within
  % an Hierarchy method and push them into this Hierarchy via this.SetLevel()
  % and/or to supply factories for automatically generating prolongators,
  % restrictors, and coarse level discretizations. Additionally contains
  % a V-cycle apply method.
  %
  % See the examples directory for numerous examples of usage.
  %

  % Note: Hierarchy inherit from VerboseObject to get verbose capability
  properties
    Levels_ % a cell-based array holding all levels in a multigrid hierarchy

    MaxCoarseSize_   = 50;        % max size of coarsest level before stopping aggregation (default: 50)
  end

  properties %(Constant)
    % default parameter
    defaultNumDesiredLevels_ = 2
  end

  methods
    function [this] = Hierarchy(arg,varargin)
      % Copy constructor
      if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
      %
      % Creates an empty multigrid Hierarchy
      this.Levels_ = {};

      % If a Matlab matrix, an operator or a level is passed
      % in arguments, fill-in the first level of the Hierarchy.
      if varexist('arg')
        if issparse(arg)
          if nargin == 1
            arg = Operator(arg);
          else
            arg = Operator(arg, varargin{:});
          end
        end
        if isa(arg,'Operator')
          FineLevel = Level();
          FineLevel.Set('A', arg); % default = SA
          arg = FineLevel;
        end
        if isa(arg,'Level')
          this.SetLevel(arg,1);
        else
          error('Hierarchy constructor: unknown arguments');
        end
      end
    end

    function SetLevel(this, level, levelId)
      % Assign a set of operators to a level within a multigrid Hierarchy
      %
      % The basic idea is that users create some of their own data and
      % operators within a multigrid Hierarchy as well as pass in
      % factories for creating other operators. The factories essentially
      % finish the specification of the multigrid Hierarchy by generating
      % all remaining operators on all levels (using the information provided
      % within user defined levels). In the normal situation
      % users only build finest grid data such as the matrix. However,
      % in a geometric method, users might provide information on coarse
      % levels as well.
      %
      % Level class describes basic functionality that a
      % level must provide. Certain multigrid methods will require
      % additional information (e.g. a Null space in smoothed aggregation).
      % See also:
      %       Level
      this.Levels_{levelId} = level;
      this.Levels_{levelId}.SetLevelId(levelId);
    end

    function [level] = GetLevel(this, levelId)
      % Retrieve operators associated with a level within a multigrid Hierarchy
      level = this.Levels_{levelId};
    end

    function Print(this)
      % Print a multigrid Hierarchy
      for i=1:length(this.Levels_)
        fprintf('******************************\n');
        fprintf('******Level %d****************\n',...
          this.Levels_{i}.GetLevelId());
        fprintf('******************************\n');
        this.Levels_{i}.Print();
      end
    end

    function SetMaxCoarseSize(this, MaxCoarseSize)
      % Determines size of coarse level when further coarsening should stop
      this.MaxCoarseSize_ = MaxCoarseSize;
    end

    function status=FullPopulate(this,Pfact, Rfact, Acfact, Smfact, startLevel,numDesiredLevels)
      % Invoke a set of factories to populate (construct prolongation,
      % restriction, coarse level discretizations, and smoothers in this
      % order) a multigrid Hierarchy starting with information on 'startLevel'
      % and continuing for at most 'numDesiredLevels'.
      %
      % Note: Empty factories are simply skipped.
      %
      % For example factories, see also Factories/Contents

      if ~varexist('startLevel'), startLevel = 1; end
      if ~varexist('numDesiredLevels'), numDesiredLevels = this.defaultNumDesiredLevels_; end

      if ~varexist('Pfact'),
          Pfact = SaPFactory();
      end
      if ~varexist('Rfact'),
          Rfact = TransPFactory(Pfact);
      end

      if ~varexist('Acfact'), Acfact = RAPFactory(); end
      if ~varexist('Smfact'),
        Smfact = SmootherFactory(Smoother('GaussSeidel', 1, 1));
        Smfact.SetOutputLevel(this.GetOutputLevel());
      end

      status = FillHierarchy(this, Pfact, Rfact, Acfact, startLevel,numDesiredLevels);

      SetSmoothers(this, Smfact, startLevel, numDesiredLevels);
    end % FullPopulate()

    function SetCoarsestSolver(this, Smfact)

      if ~varexist('Smfact'), Smfact = []; end % By default, use a direct solver

      coarsestLevelId = this.GetNumLevel();

      %% PreSetup phase (cross factory mechanism)
      if ~isempty(Smfact)
        Smfact.SetNeeds(this.Levels_{coarsestLevelId});
      end

      %% Setup phase
      this.SetupPhase(true, coarsestLevelId, 1);

      if ~isempty(Smfact)
        [Pre,Post] = Smfact.Build(this.Levels_{coarsestLevelId}); %, MySpecs);
        if ~isempty(Pre ), this.Levels_{coarsestLevelId}.Set('PreSmoother', Pre); end
        if ~isempty(Post), this.Levels_{coarsestLevelId}.Set('PostSmoother', Post); end
      else
        % TODO: this is redundant with DirectSolveSmoother.m
        % Direct solver
        level = this.Levels_{coarsestLevelId};
        Amat = level.Get('A').GetMatrixData();
        [Lfactor,Ufactor,Pperm,Qperm] = lu(Amat);
        level.Set('Lfactor',Lfactor);
        level.Set('Ufactor',Ufactor);
        level.Set('Pperm',Pperm);
        level.Set('Qperm',Qperm);
      end

      this.SetupPhase(false, coarsestLevelId, 1);
    end

    function status = SetSmoothers(this, Smfact, startLevel, numDesiredLevels)
      % Invoke a set of factories to construct smoothers within
      % a multigrid Hierarchy starting with information on 'startLevel'
      % and continuing for at most 'numDesiredLevels'.
      %
      % Note: last level smoother will not be set here. Use SetCoarsestSolver()
      % to define a smoother for the last level. Otherwise, a direct solve is
      % assumed

      if ~varexist('startLevel'), startLevel = 1; end
      if ~varexist('numDesiredLevels'), numDesiredLevels = this.GetNumLevel() - startLevel; end

      if ~varexist('Smfact')
        Smfact = SmootherFactory(Smoother('GaussSeidel', 1, 1));
        Smfact.SetOutputLevel(this.GetOutputLevel()); %TODO: should do that for every arguments ?
      end

      lastLevel = startLevel + numDesiredLevels -1;

      % Check input arguments
      if (startLevel > this.GetNumLevel()), error('startLevel to big'); end;
      if (lastLevel >= this.GetNumLevel())
        lastLevel = this.GetNumLevel()-1;
        fprintf('Warning: Coarsest Level will have a direct solve!\n');
      end

      %% PreSetup phase (cross factory mechanism)
      if ~isempty(Smfact)
        for i = startLevel:1:lastLevel,
          Smfact.SetNeeds(this.Levels_{i});
        end
      end

      %% Setup phase
      this.SetupPhase(true, startLevel, lastLevel); %numDesiredLevels;
      for i = startLevel:1:lastLevel,
        if ~isempty(Smfact),
          [Pre,Post] = Smfact.Build(this.Levels_{i}); %, MySpecs);
          if ~isempty(Pre ), this.Levels_{i}.Set('PreSmoother',  Pre ); end
          if ~isempty(Post), this.Levels_{i}.Set('PostSmoother', Post); end
        end
      end
      this.SetupPhase(false, startLevel, lastLevel); % numDesiredLevels;
      % Save stuff for outside use
      status.startLevel = startLevel;
      status.numDesiredLevels = numDesiredLevels;
      status.lastLevel = lastLevel;

      % after smoothers are set, we suppose the setup phase to be complete
      % check for memory leaks
      for i=startLevel:lastLevel % startLevel+numDesiredLevels-1
       this.Levels_{i}.TestMemoryLeak();
      end

    end % SetSmoothers()

    function status = FillHierarchy(this, Pfact, Rfact, Acfact, startLevel, numDesiredLevels)
      % Invoke a set of factories to populate (construct prolongation,
      % restriction, and coarse level discretizations in this
      % order) a multigrid Hierarchy starting with information on 'startLevel'
      % and continuing for at most 'numDesiredLevels'.
      %
      % For example factories, see also Factories/Contents
      %
      if ~varexist('startLevel'), startLevel = 1; end
      if ~varexist('numDesiredLevels'), numDesiredLevels = this.defaultNumDesiredLevels_; end

      if this.GetOutputLevel() > 0
        fprintf('Hierarchy: start level     = %d\n',startLevel);
        fprintf('Hierarchy: maximum #levels = %d\n',numDesiredLevels);
      end

      if varexist('Rfact') && varexist('Pfact') && isempty(Rfact)
          if ~isa(Pfact,'PRFactory')
              warning('Rfact is empty, but Pfact is not an PRFactory');
          end
          PRfact = Pfact; % note: max coarse size of PRfact won't be overwritten!
      elseif varexist('Rfact') && varexist('Pfact') && ~isempty(Rfact)
        PRfact = GenericPRFactory(Pfact,Rfact);
        PRfact.SetMaxCoarseSize(this.MaxCoarseSize_);
      elseif varexist('Pfact') && ~varexist('Rfact')
        Rfact = TransPFactory(Pfact); %SaPFactory();
        PRfact = GenericPRFactory(Pfact,Rfact);
        PRfact.SetMaxCoarseSize(this.MaxCoarseSize_);
      else
        PRfact = GenericPRFactory(SaPFactory);
        PRfact.SetMaxCoarseSize(this.MaxCoarseSize_);
      end
      if ~varexist('Acfact'), Acfact = RAPFactory(); end

      % Check method arguments
      if (startLevel > length(this.Levels_)) || (isempty(this.Levels_{startLevel})), error('Level %d does not exist. Must be created before calling FillHierarchy()', startLevel); end
      if startLevel ~= this.Levels_{startLevel}.GetLevelId(), error('Level %d return GetLevelId() == %d', startLevel, this.Levels_{startLevel}.GetLevelId()); end %TODO: useless ? Or check every levels ?

      % Create the level hierarchy
      % Note: this.Levels_{startLevel} already exist and we create this.Levels_{startLevel+1}, +2... The last created level is: this.Levels_{(numDesiredLevels-1)+1}
      for i=startLevel:numDesiredLevels-1
        if (i+1 > length(this.Levels_)) || (isempty(this.Levels_{i+1}))
          this.Levels_{i+1} = this.Levels_{i}.BuildMe();
          this.Levels_{i+1}.SetLevelId(i+1);
        end
      end

      %% PreSetup phase (cross factory mechanism)
      for i=startLevel:numDesiredLevels-1
        if ~isempty(PRfact), PRfact.SetNeeds(this.Levels_{i}, this.Levels_{i+1}); end
        if ~isempty(Acfact), Acfact.SetNeeds(this.Levels_{i}, this.Levels_{i+1}); end
      end

      %% Setup phase
      this.SetupPhase(true, startLevel, numDesiredLevels);

      for i=startLevel:numDesiredLevels-1
        if ~isempty(PRfact)
          flag = PRfact.Build(this.Levels_{i},this.Levels_{i+1}); %, MySpecs);
          if ~flag, this.Levels_ = {this.Levels_{1:i}}; break; end % note: {} seems mandatory
        end

        if ~isempty(Acfact)
          flag=Acfact.Build(this.Levels_{i},this.Levels_{i+1}); %,MySpecs);
          if ~flag, this.Levels_ = {this.Levels_{1:i}}; break; end % note: {} seems mandatory
        end
      end

      % End of setup phase
      numLevel = this.GetNumLevel(); this.SetupPhase(false, startLevel, numLevel);

      if this.GetOutputLevel() > 0
        fprintf('Hierarchy: actual #levels  = %d\n', numLevel);
      end

      % TMP DEBUG Test (should be done after smoother in fact, at
      % the very end)
      %for i=startLevel:startLevel+numLevel-1
      %  this.Levels_{i}.TestMemoryLeak();
      %end

      status.lastLevel  = this.GetNumLevel();
    end

    function status = GetStatistics(this)

      %% Statistics
      status.lastLevel  = this.GetNumLevel();

      % Compute FineNnz (TODO: if this.Levels_{i}.IsAvailable('A') ?)
      startA = this.Levels_{1}.Get('A').GetMatrixData();
      if issparse(startA), status.FineNnz = nnz(startA); else status.FineNnz  = -1; end

      % Compute TotalNnz
      status.TotalNnz = 0;
      for i=1:this.GetNumLevel()
        if this.Levels_{i}.IsAvailable('A')
          A = this.Levels_{i}.Get('A').GetMatrixData();
          if issparse(A)
            status.TotalNnz = status.TotalNnz + nnz(A);
          else
            status.TotalNnz = -1;
            break;
          end
        end
      end

      % Compute Operator Complexity
      if status.FineNnz ~= -1 && status.TotalNnz ~= -1, status.OperatorComplexity = status.TotalNnz/status.FineNnz; else status.OperatorComplexity = -1; end

    end

    % sometime, we need to know how many levels there are.
    function numlevel = GetNumLevel(this)
      numlevel = length(this.Levels_);
    end

    function [sol] = Iterate(this, rhs, Nits, sol, InitGuessStatus, Cycle, levelId)
       % Perform one iteration of a  multigrid cycle
       % starting from level levelId (if supplied, otherwise start from level 1).
       %
       % This is a recursive method
       %
       % Note: Optionally, InitGuessStatus can be passed in. If 0, it
       %       indicates that sol is initially zero and so this can save
       %       some work in the computations.

       mue_include
       if ~varexist('levelId'),    levelId = 1; end
       if ~varexist('Cycle'),      Cycle = 1; end
       if ~varexist('InitGuessStatus'), InitGuessStatus = NOTALLZEROS; end
       if ~varexist('sol'), sol = zeros(size(rhs,1),1); InitGuessStatus = ALLZEROS; end %TODO:multiple_rhs

%        % DEBUG
%        % check this level for memory leaks.
%        % not really necessary, since it's already done in SetSmoothers
%        this.Levels_{levelId}.TestMemoryLeak();
%
       for i=1:Nits
          Fine = this.Levels_{levelId};

          if length(this.Levels_) ~= levelId
             % Recursive multigrid: Pre Smoothing, Coarse grid correction and Post-Smoothing

             Coarse       = this.Levels_{levelId+1};

             % Pre-Smoothing
             if Fine.IsAvailable('PreSmoother')
                [sol,InitGuessStatus]= Fine.Get('PreSmoother').Apply(sol,rhs,InitGuessStatus);
             end

             % Residual
             r   = rhs - Fine.Get('A') * sol;

             if Fine.GetLevelId()==1 && this.GetOutputLevel() > 0,
                fprintf('%3d: ||r||=%e\n',i,norm(r));
             end

             % Restrict residual and solution to the coarse grid
             rhat = Coarse.Get('R') * r;

             [n,m] = size(rhat); % This is a FIX: See commit e95579f56f58 for more info. To be investigated...
             if (isa(rhat, class('MultiVector'))), csol = MultiVector(n,m); else csol = zeros(n,m); end

             % Recursive call
             csol = this.Iterate(rhat,1,csol, ALLZEROS, Cycle, levelId+1);
             if Cycle>1
                csol = this.Iterate(rhat,1,csol, NOTALLZEROS, Cycle,levelId+1);
             end

             % New fine solution
             sol = sol + Coarse.Get('P') * csol;
             InitGuessStatus = NOTALLZEROS;

             % Post-Smoothing
             if Fine.IsAvailable('PostSmoother')
                [sol,InitGuessStatus]= Fine.Get('PostSmoother').Apply(sol,rhs,InitGuessStatus);
             end

          else % End of the recursion

             % Enable a faster coarse solver - TODO: should be done
             % somewhere else
             if ~Fine.IsAvailable('PreSmoother') && ~Fine.IsAvailable('PostSmoother')
                 this.SetCoarsestSolver(); % version 1
                 this.SetCoarsestSolver(SmootherFactory(DirectSolveSmoother,[])); % version 2
             end

             if Fine.IsAvailable('PreSmoother')
                [sol,InitGuessStatus] = Fine.Get('PreSmoother').Apply(sol,rhs,InitGuessStatus);
             end
             if Fine.IsAvailable('PostSmoother')
                [sol,InitGuessStatus] = Fine.Get('PostSmoother').Apply(sol,rhs,InitGuessStatus);
             end

             if ~Fine.IsAvailable('PreSmoother') && ~Fine.IsAvailable('PostSmoother')

                if Fine.IsAvailable('Lfactor') && Fine.IsAvailable('Ufactor') && Fine.IsAvailable('Pperm') && Fine.IsAvailable('Qperm')
                   sol = Fine.Get('Qperm')*(Fine.Get('Ufactor')\(Fine.Get('Lfactor')\(Fine.Get('Pperm')*rhs)));

                else
                   % TODO: eventually, all of the following should go in a default Direct Solve
                   % Smoother class...

                   if 1
                       sol = Fine.Get('A') \ rhs; % slow!! direct factorization done at each iteration
                   else
                       % The following code had been introduced for singular matrices.
                       % TODO: We might want to make this a requestable option for singular matrices
                       % Note this version doesn't work for size(oo,2) > 1.

                       temp = 1;
                       if Fine.IsAvailable('NullSpace')
                           % TODO: The qr factorization and 'temp' should not be recomputed at each iteration.
                           oo  = Fine.Get('NullSpace'); [qq,rr] = qr(oo,0); oo = qq;
                           A =  Fine.Get('A').GetMatrixData;
                           Anorm = norm(A,'fro');
                           if size(A,2) == size(oo,1), temp = (A * oo)/Anorm; end;
                       end

                       if norm(temp,'fro') > 1e-12,
                           sol = Fine.Get('A') \ rhs;
                       else
                           aaa = Fine.Get('A').GetMatrixData; nnn = size(aaa,1);
                           sol = zeros(nnn,1);
                           minor = (1:nnn-size(oo,2));  %
                           sol(minor) = aaa(minor,minor)\rhs(minor);
                           sol = sol - (sol'*oo)*oo;
                       end
                   end % if 1

                end
             end
             if Fine.GetLevelId() == 1 && this.GetOutputLevel() > 0,
                r = rhs - Fine.Get('A')*sol;
                fprintf('%3d: ||r||=%e\n',i,norm(r));
             end

          end
       end
    end

    function TestMemoryLeak(this)
      for i=1:this.GetNumLevel()
        this.Levels_{i}.TestMemoryLeak();
      end
    end

  end

  methods (Access = private)
    function SetupPhase(this, TorF, startLevel, numDesiredLevels)
      for i=startLevel:numDesiredLevels
        this.Levels_{i}.SetupPhase(TorF);
      end
    end
  end

  methods (Access = protected)

    function Copy_(this, src, mc)
      %COPY_
      %
      %   SYNTAX   obj.Copy_(src, mc);
      %
      %     src - Object to copy
      %     mc  - MATLAB Metaclass
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);
    end

  end % methods

end
