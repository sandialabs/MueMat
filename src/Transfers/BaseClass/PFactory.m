%% PFactory
% interface class for prolongation operator
%%
classdef PFactory < VerboseObject
    % This factory provides an interface for a concrete implementation of a
    % prolongation operator (e.g. tentative prolongator)
    % For a concrete implementation the user has to overwrite the virtual
    % Build method

   properties (Access = protected)
        prolongation_mode_ = true; % default: PgPFactory is in prolongation mode
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% public functions
   methods
       function [this] = PFactory(arg)
       %PFACTORY Constructor

           % copy constructor
           if nargin == 1 && isa(arg,class(this)), this.Copy_(arg,[]); return; end;
       end

       function [ToF] = ProlongationMode(this, ToF)
          if varexist('ToF'),
              ToFold = this.prolongation_mode_;
              this.prolongation_mode_ = ToF;
              ToF = ToFold;
          else ToF = this.prolongation_mode_; end
       end

       function [flag] = CheckForReUsableP(this, FineLevel, CoarseLevel, Pfact)
           flag = false;

           if this.ProlongationMode() == true
               if CoarseLevel.IsAvailable('P', Pfact),
                    flag = true;
               end
           else
               if CoarseLevel.IsAvailable('R', Pfact),
                    flag = true;
               end
           end
       end
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% abstract functions
   methods (Abstract)
     flag = Build(this, FineLevel, CoarseLevel, Specs);

     SetNeeds(this, FineLevel, CoarseLevel);
     % To place request on the data of the level (increment counter
     % by calling Level.Request())

     [ToF] = SupportsRestrictionMode(this);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% protected functions
   methods (Access = protected)

       function Copy_(this,src,mc)
           % COPY_
           % syntax obj.Copy_(src,mc);
           % src: object to copy
           % mc: MATLAB Metaclass
           [cmd,data,mc] = this.CopyCmd_(src,mc);
           eval(cmd);
       end
   end

end
