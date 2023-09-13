% A factory which creates aggregates in the structure AggInfo
%
classdef SillyAggFact < AggregationFactory
   properties (Access = private)
      Needs_    % indicates cross-factory needs
   end
   methods
      function [this] = SillyAggFact()
      % create a silly aggregation factory that requests a graph to be saved.
      % this class is only meant to test some reuse options.
         this = this@AggregationFactory();
         Need.SaveAggregates='all';
         this.AddNeeds(Need);
         fprintf('SillyAggFact: Need Graph\n');
      end
   end
end
