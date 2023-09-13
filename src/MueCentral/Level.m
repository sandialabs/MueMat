classdef Level <  Needs

  % -- old comments, to be updated ...
  %   % Base class for data associated with a single level of a multigrid hierarchy
  %   % All data buckets must provide methods to set and get a prolongator,
  %   % a restrictor, a discretization operator, and a smoother. In addition,
  %   % there must be a means of setting and getting an id number associated
  %   % with the level.
  %   %
  %   % Data associated with one level of a smoothed aggregation method
  %   % In addition to standard AMG operators, this data bucket also provides
  %   % a near null space associated with the discretization operator.
  % --

  % Commonly used data:
  % A             % discretization operator
  % Afiltered     % modified version of A  with some entries dropped
  %               % and perhaps the diagonal modified to maintain rowsums
  % R             % restriction operator
  % P             % prolongator operator
  % PreSmoother   % smoother operator
  % PostSmoother  % smoother operator
  % LevelId       % id number associated with level
  %
  % Smoothed Aggregation data:
  % NullSpace     % Near null space of discretization operator
  % AggInfo       % Aggregation information which can be optionally saved
  % Ptxnt         % tentative prolongator which can be optionally saved
  % Rtent         % for unsymmetric problemes
  % Graph         % saving a graph for SA
  %
  % Less common data:
  % AuxMatrix     % can be used instead of A  to steer coarsening/sparsity pattern algorithms
  % AuxMatrixFunc % function used to regenerate AuxMatrix from A & coords instead of via RAP
  % AuxMatP       % an alternative prolongator for projecting AuxMatrix
  % AuxMatR       % an alternative restrictor for projecting AuxMatrix
  % xcoords       % coordinates can be used to steer coarsening/sparsity or for visualization
  % ycoords       % coordinates can be used to steer coarsening/sparsity or for visualization
  % zcoords       % coordinates can be used to steer coarsening/sparsity or for visualization
  % coordP        % an alternative prolongator for projecting coordinates
  % coordR        % an alternative restrictor for projecting coordinates

  properties (Access = private)
    LevelId_      % id number associated with level
  end

  methods %(Static = true) not static anymore (wasn't use as a static method by the way)
    function [newLevel] = BuildMe(this)
      % This method is used by the class Hierarchy to create additional MG
      % levels when such levels are needed and are not directly created by users.
      % In Hierarchy.m, the last level i of the hierarchy is copied
      % to create a new level i+1.
      %
      % If you use your own Level class (a specialization of this
      % one), you should specialized this method. Then,
      % Hierarchy will use automatically the right kind of Level in
      % your hierarchy.

      % The keep/keepAll status of the current level i is
      % also copied to the level i+1.
      %
      % Why? With this default behavior, you can use easily use the
      % 'reuse' capability of MueMat by providing only the FineLevel:
      %
      % In the following example:
      %  FineLevel = Level();
      %  FineLevel.Keep('P'); % disable the auto desallocation of
      %                       % 'P' to be able to reuse it.
      %  MgHierarchy.SetLevel(FineLevel,1);
      %
      % the 'Keep' option of 'P' will be propagate to each level of
      % your hierarchy.
      %
      % To disable this default behavior, define every level
      % of your hierarchy manually or overload this class.

      % Create an instance.
      newLevel = Level();

      % Copy keepAll status
      newLevel.KeepAll(this.IsKeptAll());

      % Copy independant keep status
      for k=keys(this.countTable_)
        ename = k{1};

        for l = this.countTable_.handles(ename)
            ehandle = l{1};
            if this.IsKept(ename, ehandle)
              newLevel.Keep(ename, ehandle);
            end
        end
      end

    end
  end

  methods
    function [this] = Level(arg)
      % Copy constructor
      if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
      %

      this.Keep('A');
      this.Keep('P');
      this.Keep('R');

      this.Keep('PreSmoother');
      this.Keep('PostSmoother');
    end

    function SetLevelId(this,i)
      % Assigns an id to this bucket
      this.LevelId_ = i;
    end

    function [z] = GetLevelId(this)
      % Retrieves id associated with this bucket
      z = this.LevelId_;
    end

    % Print()
    function Print(this)
      fprintf('LevelId: %d\n', this.LevelId_);
      Print@Needs(this);
    end

  end % methods

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

end %class
