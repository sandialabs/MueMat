classdef NeedsObject < VerboseObject
% Base class that provides Needs list.
% Needs are simply data that an object requires but that might be created
% elsewhere.  The object registers that it requires the data via the Needs
% list.  For example, a restriction factory that transposes the tentative
% prolongator 'Needs' the prolongator factory to save this.
   properties (Access = private)
      Needs_          % Structure holding needs (see CrossFactory)
   end

   methods
      function [this] = NeedsObject(arg)
         % Copy constructor
        if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
        %

         this.Needs_ = [];
      end
      function AddNeeds(this,Needs)
      % add a new need to the factory's list of cross-factory specifications
          this.Needs_ =CrossFactory.MergeNeeds(Needs,this.Needs_);
      end
      function [z] = GetNeeds(this)
      % retrieve all needs of the factory.
          z = this.Needs_;
      end

  end %public methods

  methods (Access = protected)

    function Copy_(this, src, mc)
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);
    end

  end % methods

end % classdef
