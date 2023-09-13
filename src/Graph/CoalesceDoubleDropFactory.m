classdef CoalesceDoubleDropFactory < CoalesceDropFactory
% A variant of the CoalesceDropFactory that does two passes
% of dropping.  The first dropping is done on the amalgamated or
% AuxMatrix.  The second is done on the original matrix as filtered
% by dropping on the amalgamated/aux one.
%

% TODO: this design is very confusing. One factory should be used for each
% matrix that has to be coalesce and/or dropped. CoalesceDoubleDropFactory
% should not call Build@CoalesceDropFactory(this,CurrentLevel);

    properties (Access = private)
     DoublePreDropFunc_  = [];
     DoublePreDropData_  = [];    % User supplied data passed to DoublePreDropFunc_()
     DoublePostDropFunc_ = [];
     DoublePostDropData_ = [];    % User supplied data passed to PreDropFunc_()
   end
   methods

      function [this] = CoalesceDoubleDropFactory(arg)
         % Copy constructor
         if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
         %
      end

      function SetDoublePreDropSpecifications(this,DoublePreDropFunc,DoublePreDropData)
      % Assign predrop function used to decide matrix elements ignored before amalgamation occurs
         this.DoublePreDropFunc_ = DoublePreDropFunc;
         if varexist('DoublePreDropData'), this.DoublePreDropData_ = DoublePreDropData; end
      end
      function SetDoublePostDropSpecifications(this,DoublePostDropFunc,DoublePostDropData)
      % Assign postdrop function used to decide matrix elements ignored after amalgamation occurs
         this.DoublePostDropFunc_ = DoublePostDropFunc;
         if varexist('DoublePostDropData'), this.DoublePostDropData_     = DoublePostDropData; end
      end

      function [DoublePreDropFunc, DoublePreDropData] = GetDoublePreDropSpecifications(this)
         DoublePreDropFunc = this.DoublePreDropFunc_;
         DoublePreDropData = this.DoublePreDropData_;
      end
      function [DoublePostDropFunc, DoublePostDropData] = GetDoublePostDropSpecifications(this)
         DoublePostDropFunc = this.DoublePostDropFunc_;
         DoublePostDropData = this.DoublePostDropData_;
      end
      %function SetNeeds(this, CurrentLevel)
      %    SetNeeds@CoalesceDropFactory(CurrentLevel);
      %end

      function flag = Build(this,CurrentLevel,Specs)

      % Secondary PreDrop
      % This will have to force regeneration of the AuxMatrix if
      % we're using one.

      % Primary CoalesceDrop
      Build@CoalesceDropFactory(this,CurrentLevel);

      % Get the name of the A matrix
      % NTS: might want to make it user controllable.
      A=CurrentLevel.Get('A');

      % Build an amalgamated A
      this.SetBinary();
      [AmalgFunc, AmalgData] = this.GetAmalgSpecifications();
      Aamalg=this.AmalgamateAndOrDrop(A,CurrentLevel,[],[],AmalgFunc,AmalgData,[],[],'NotBinary');

      % Filter the amalgamated A
      Af=CurrentLevel.Get(this.GetAFilteredName()).GetMatrixData().*Aamalg.GetMatrixData;
      Aamalg.SetMatrixData(Af);
      CurrentLevel.Set('SecondaryA',Aamalg);

      % Secondary PostDrop
      [Graph,Afiltered]=this.AmalgamateAndOrDrop(Aamalg,CurrentLevel,[],[],AmalgFunc,AmalgData,this.DoublePostDropFunc_,this.DoublePostDropData_,'Binary');

      CurrentLevel.Set('Graph', Graph);
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

end %classdef CoalesceDropFactory

