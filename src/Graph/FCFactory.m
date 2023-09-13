% A factory which creates aggregates in the structure AggInfo
%
classdef FCFactory < VerboseObject
   properties (Access = private)
      Algorithm_
      TargetSize_
      PointsPerDim_
   end
   methods
      function [this] = FCFactory(obj)
         % Copy constructor
         if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
         %

         % create an FC factory
         this.Algorithm_ = 'graph';
         this.TargetSize_ = 14;
      end
      function [FCPts] = Build(this, Graph, Specs)
      % take a matrix graph and build a set of F and C points
      %
      % On output:
      %     FCPts.Fpts              list of Fpoints
      %
      %     FCPts.Cpts              list of Cpoints
      %
      % where  Union(FCPts.Fpts,FCPts.Cpts) gives all vertices and
      %        Intersection(FCPts.Fpts,FCPts.Cpts) is empty.
      %
        this.TempOutputLevel(Specs);
        [AggregateId, roots] = AggregationFactory.GraphAggregation(Graph.MatrixData);
        FCPts.Cpts = roots;
        temp = ones(Graph.RowMap.NNodes(),1);
        temp(roots) = 0;
        FCPts.Fpts = find(temp);
        this.RestoreOutputLevel();
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
