classdef VerboseObject < CopiableHandle
% Base class to control output.  All other classes derive from this one.
   properties (Access = private)
      OutputLevel_    % indicates level of output
      PriorOutputLevel_ = [];
   end
   methods
      function [this] = VerboseObject(arg)
         % Copy constructor
        if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
        %

         this.OutputLevel_ = 0;
      end
      function SetOutputLevel(this,OutputLevel)
      % Set the output level
        this.OutputLevel_ = OutputLevel;
      end
      function [OutputLevel] = GetOutputLevel(this)
      % get the output level
        OutputLevel = this.OutputLevel_;
      end

      function TempOutputLevel(this,Specs)
      % used to temporarily set Output Level during build() invocation if set via a cross-factory specification
        this.PriorOutputLevel_ = this.OutputLevel_;
        if ~isempty(Specs)
           temp = Specs.GetSpecifications();
           if isfield(temp,'OutputLevel'),
              this.OutputLevel_ = max(this.PriorOutputLevel_,temp.OutputLevel);
           end
        end
      end

      function RestoreOutputLevel(this)
      % used to restore Output Level to state before invocation of TempOutputLevel()
         if isempty(this.PriorOutputLevel_)
            fprintf('VerboseObject::TempOutputLevel not invoked prior to VerboseObject::RestoreOutputLevel\n');
            keyboard;
         end
         this.OutputLevel_ = this.PriorOutputLevel_;
         this.PriorOutputLevel_ = [];
      end

  end

  methods (Access = protected)

    function Copy_(this, src, mc)
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);
    end

  end % methods

end % classdef
