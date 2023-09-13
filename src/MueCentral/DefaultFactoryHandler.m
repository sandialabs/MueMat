classdef DefaultFactoryHandler < VerboseObject
   % class that provides default factories within Needs class

   properties (Access = private)
        factoryTable_;  % hashtable: variable name -> factory class
   end

   methods
        function [this] = DefaultFactoryHandler(varargin)
            if nargin == 1 && isa(varargin, class(this)), this.Copy_(varargin,[]); return; end

            this.factoryTable_ = Hashtable();

            % add some default factories
            this.SetDefaultFactory('P',TentativePFactory());
%            this.SetDefaultFactory('Aggregates',AggregationFactory()); % with default CoalesceDropFactory
        end

        function SetDefaultFactory(this,varname,ehandle)
           if ~ischar(varname)
              error('varname must be a string\n');
           end
           if ~ismethod(ehandle,'Build')
              error('ehandle is not a factory class\n');
           end
           this.factoryTable_.Set(varname,ehandle);
        end

        function [ehandle] = GetDefaultFactory(this,varname)
           if ~ischar(varname)
              error('varname must be a string\n');
           end
           if ~this.factoryTable_.isKey(varname)
              error('no default factory for %s',varname);
           end

           ehandle = this.factoryTable_.Get(varname);
        end
   end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%
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