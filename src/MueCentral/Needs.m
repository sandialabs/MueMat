classdef Needs < VerboseObject
    %
    % Class that allows cross-factory communication of data needs.
    %
    % Maintains a list of 'Needs' for a given Level. For example, a restriction factory that
    % transposes the tentative prolongator 'Needs' the prolongator factory to save this.
    %
    % TODO: update the description of this file
    % Use a reference counts to deallocate objects which are no longer referenced

    properties (Access = protected)
        countTable_; % Stores number of outstanding requests for a need.
        dataTable_;  % Stores data associated with a need.

        keepAll_ = 1 % Keep all the data. No desallocation using reference counting (for reuse or debug)
        debug_   = 0 %

        setupPhase_ = 0 % 0 Set() always put data (used by user before presetup phase)
        % 1 Set() put data only iff countTable(key) exist and is ~=0

        defaultFactoryHandler_ = DefaultFactoryHandler();
    end % properties

    methods

        function [this] = Needs(varargin)
            if nargin == 1 && isa(varargin, class(this)), this.Copy_(varargin,[]); return; end

            this.countTable_ = ExtendedHashtable();
            this.dataTable_  = ExtendedHashtable();
        end

        %% public helper functions

        function SetDefaultFactoryHandler(this, handler)
           this.defaultFactoryHandler_ = handler;
        end

        function [handler] = GetDefaultFactoryHandler(this)
           handler = this.defaultFactoryHandler_;
        end

        % helper functions, that translates input ename and ehandle_in into
        % a valid factory handle
        % ehandle_in can be
        %      'default' (string): call this.defaultFactoryHandler_ for
        %                          to obtain default factory for variable
        %                          ename, if available
        %       0        (int)   : no generating factory available, ename is
        %                          just a plain variable -> do nothing
        %       handle (Factory) : do nothing, ehandle_in is already a
        %                          (valid) factory
        function [ehandle_out] = InterpretHandle(this, ename, ehandle_in)
           if ischar(ehandle_in)
                if strcmp(ehandle_in,'default')
                    ehandle_out = this.defaultFactoryHandler_.GetDefaultFactory(ename);
                else
                    fprintf('ehandle %s not recognized\n',ehandle_in);
                    error('invalid ehandle');
                end
           else
                ehandle_out = ehandle_in;
           end
        end

        %% Manage memory allocation / desallocation using a reference counter mecanism.

        % Indicate that an object is needed. This increments the storage counter.
        % This method is used during the pre-setup phase.
        function Request(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);

            % Following line for debugging:
            %if strcmp(ename,'Aggregates'), fprintf('Needs request: %s\n', ename); dbstack; end

            % If it the first request, create a new key in the hashtable
            if ~this.countTable_.isKey(ename, ehandle), this.countTable_.Set(ename, 0, ehandle); end;

            if this.debug_, this.countTable_.Set(ename,-1,ehandle); end; % disable reference counting

            % Increment the counter
            this.IncrementCounter(ename, ehandle);
        end % Request()

        % Decrement the storage counter associated with ename.
        % This method is used during the setup phase.
        % This method use both countTable_ and dataTable_
        function Release(this,ename,ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);

            % Test: we can only release data if key exists on countTable...
            if ~this.countTable_.isKey(ename,ehandle) || ~this.dataTable_.isKey(ename,ehandle), fprintf('Needs.Release(): %s not found\n',ename); keyboard; end

            % Decrement the reference counter
            this.DecrementCounter(ename, ehandle);

            % Desallocation if counter reachs zero.
            if this.countTable_.Get(ename,ehandle) == 0 && this.keepAll_ == false
                this.countTable_.Remove(ename,ehandle);
                this.dataTable_.Remove(ename,ehandle);
            end
        end % Release()

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        %% Set() and Get() methods.
        % The following methods doesn't modify the reference counter and work
        % only with dataTable_.

        % Store data. This does not modify the storage counter.
        function Set(this, ename, entry, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);

            %if this.IsAvailable(ename), warning('MueLu:LevelSetOverwrite','Set(): data for key ''%s'' already present in the Level (will be overwritten)', ename); end

            if this.debug_, this.countTable_.Set(ename,-1,ehandle); end; % disable reference counting

            if ~this.setupPhase_
                % Always store data in this mode
                if ~this.countTable_.isKey(ename,ehandle), this.countTable_.Set(ename,0,ehandle); end
                this.dataTable_.Set(ename,entry,ehandle);
            else
                % We store data iff data are needed by another factory, ie: counter~=0
                if (this.countTable_.isKey(ename,ehandle) && this.countTable_.Get(ename,ehandle) ~= 0) || this.keepAll_ == true
                    if ~this.countTable_.isKey(ename,ehandle), this.countTable_.Set(ename,0,ehandle); end % it can happen for keepAll == true
                    this.dataTable_.Set(ename,entry,ehandle);
                end
            end

        end % Set()

        % Get data. This does not decrement the storage counter.
        function data = Get(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);

            if ~this.dataTable_.isKey(ename,ehandle),
                % do not just throw an error, but first check if - even though a
                % valid ehandle is missing, to find a unique answer to the
                % request for variable ename
                fprintf('Needs.Get(): %s not found. Did you call %s.Build before?\n',ename, ehandle);
                % first of all we know, that ename is a valid key for dataTable_
                % then check the number of keys of generating factories for
                % ename
                keylist = this.dataTable_.handles(ename);
                if length(keylist) > 1
                    % bad: there are more than one requests for ename with differing
                    % generating factories, that is we can not find a unique
                    % answer to the request for ename without a valid ehandle
                    error('There are more than one request for ''%s'' with differing ehandles (generating factories). You must also provide an instance of the factory that created this data.',ename);
                else
                    warning('return unique answer for %s! This code is not safe and may cause problems! Avoid requesting variables without valid ehandle!\n',ename);
                    ehandle = keylist{1};
                end

                % we cannot build data within level class, since it could
                % be, that we need Fine and CoarseLevel :-(
                % check if ehandle is a valid factory
                %if ~ismethod(ehandle,'Build'), error('ehandle %s is not a valid factory', ehandle); end;
            end
            data = this.dataTable_.Get(ename,ehandle);
        end % Get()

        % Test whether a need's value has been saved.
        function [status] = IsAvailable(this,ename,ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            status = this.dataTable_.isKey(ename,ehandle);
        end % IsAvailable()

        % Test whether a need's value is Requested.
        function [status] = IsRequested(this,ename,ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            status = this.countTable_.isKey(ename,ehandle);
        end % IsRequested()
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        %% Utilities

        % Switch Setup Phase
        function SetupPhase(this, TorF)
            this.setupPhase_ = TorF;

            if this.setupPhase_ == true
                if ~this.keepAll_
                    % Desallocation of data that have not been requested
                    % JG TODO: problem: if dta needed by smoothers!!!!
                    for k=this.dataTable_.keys();
                        ename = k{1};

                        for l=this.dataTable_.handles(ename);
                            ehandle = l{1};
                            if ~this.countTable_.isKey(ename,ehandle) || this.countTable_.Get(ename,ehandle) == 0
                                warning('MueMat:NeedsSetupPhaseDesalloc','''%s'' is present in the Level structure but have not been requested (=> removed)', ename);
                                this.dataTable_.Remove(ename,ehandle);
                                this.countTable_.Remove(ename,ehandle);
                            end
                        end
                    end
                end
            end
        end

        % Print()
        function Print(this)

            % Basic output :
            %  fprintf('Counters:\n'); %disp(this.countTable_);
            %  disp(keys(this.countTable_));
            %  disp(values(this.countTable_));
            %  fprintf('Data:\n');     %disp(this.dataTable_);
            %  disp(keys(this.dataTable_));

            fprintf('Options:\n');
            fprintf('    keepAll    = %d\n', this.keepAll_);
            fprintf('    debug      = %d\n', this.debug_);
            fprintf('    setupPhase = %d\n', this.setupPhase_);
            fprintf('\n');

            fprintf('Data: \n');
            fprintf('    Key []                     Counter            Value\n\n');

                  for k=this.countTable_.keys()
                     ename = k{1};

                     for l=this.countTable_.handles(ename);
                         ehandle = l{1};

                         if this.dataTable_.isKey(ename, ehandle)
                            % dataStr = 'X';
                            sizeStr = sprintf('%dx', size(this.dataTable_.Get(ename,ehandle))); sizeStr = sizeStr(1:end-1);
                            dataStr = sprintf('%s %s', class(this.dataTable_.Get(ename,ehandle)), sizeStr);
                         else
                            dataStr = ' ';
                         end

                         % compute spacing
                         nSpace = max(20 - length(ename) - length(class(ehandle)),2); spaceStr = repmat(' ',1,nSpace);
                         nSpace2 = max(15 - length(sprintf('%d',this.countTable_.Get(ename,ehandle))),2); spaceStr2 = repmat(' ',1,nSpace2);

                         fprintf('    ''%s'' [%s] %s [%d] %s %s\n', ename, class(ehandle), spaceStr, this.countTable_.Get(ename,ehandle), spaceStr2, dataStr);
                     end
                  end

                  fprintf('\n');



         end % Print

    function TestMemoryLeak(this)
      % This method checks that no temporary data remains on the
      % hashtable. Keeped data are not considered as temporary data.

      % Note that this method do not make any assumption about the sync of
      % countTable and dataTable and thus can be able to catch out of
      % sync bug.

      % Following lines for debug:
      % warning('MueMat:TestMemoryLeak()');
      % this.Print();

      if ~this.keepAll_
        for k=keys(this.countTable_)
          ename = k{1};
          for l=this.countTable_.handles(ename);
              ehandle = l{1};
              if this.countTable_.Get(ename,ehandle) ~= -1
                if this.dataTable_.isKey(ename,ehandle)
                  warning('MueMat:NeedsTestMemoryLeak','Key and Data ''%s'' (handle %s) still present in the data structure', ename, class(ehandle));
                else
                  warning('MueMat:NeedsTestMemoryLeak','Key          ''%s'' (handle %s) still present in the data structure', ename, class(ehandle));

                end
            end
        end
        end
      end
    end
    end % methods

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    %% Disable the reference counter mecanism
    % The couple Keep()/Delete() allows a user to specify that he want to keep
    % the data after the run (for reuse in another run, plotting or debug).
    % User is then responsible to delete the data.
    % The current implementation use a special value (-1) of the
    % counter to disable it.

    methods (Access = public)

        function Keep(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            % TODO: warning if hashtable entry already exist ?
            this.countTable_.Set(ename,-1,ehandle);
        end

        function KeepAll(this, TorF)
            if ~varexist('TorF'), TorF = true; end
            this.keepAll_ = TorF;

            % TODO: when KeepAll is turned off, remove data for counter <= 0
        end

        function TorF = IsKept(this,ename,ehandle)
            % note: give only the personal Keep status of ename. Don't take
            % into account the global keep all option. Use IsKeptAll() for that
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            if ~this.countTable_.isKey(ename,ehandle), TorF = false; return; end;
            TorF = this.countTable_.Get(ename,ehandle) == -1;
        end

        function TorF = IsKeptAll(this)
            TorF = this.keepAll_;
        end

        function Delete(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            if ~this.countTable_.isKey(ename,ehandle), fprintf('Needs.Delete(): %s not found\n',ename); keyboard; end

            if (this.countTable_.Get(ename,ehandle) ~= -1) && this.keepAll_ == false
                fprintf('this.Delete(): This method is intend to be use when the automatic garbage collector is disable (using this.Keep() or this.KeepAll()). To decrement the reference counter of ''%s'', please use this.Release()\n',ename);
                keyboard;
            end

            this.countTable_.Remove(ename,ehandle);
            if this.dataTable_.isKey(ename,ehandle), this.dataTable_.Remove(ename,ehandle); end % if the data is not present, it's a silent error.
        end

    end %methods
    methods (Access = private)

        function IncrementCounter(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            if (this.countTable_.Get(ename,ehandle) ~= -1) % handle the case where the counter is disable
                this.countTable_.Set(ename,this.countTable_.Get(ename,ehandle) + 1,ehandle);
            end
        end

        function DecrementCounter(this, ename, ehandle)
            if ~varexist('ehandle'), ehandle = 0; end;
            ehandle = this.InterpretHandle(ename,ehandle);
            if (this.countTable_.Get(ename,ehandle) ~= -1) % handle the case where the counter is disable
                cnt = this.countTable_.Get(ename,ehandle);
                this.countTable_.Set(ename,cnt-1,ehandle);
            end
        end

    end % methods

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

end % class Needs

%TODO:
% - autodetect out of sync between the two hastable
% - what happen if you turn on Keep or KeepAll and then turn them off ? =>
%   should we desallocate data with -1 (for Keep) and counter <0 (for KeepAll) ?
