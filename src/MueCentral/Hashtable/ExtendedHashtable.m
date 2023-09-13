classdef ExtendedHashtable < CopiableHandle
% ExtendedHashtable
%
% Class provides MATLAB work around for a hash table in a hashtable with
% MATLAB handles as keys

  properties (Access = private)
    dataTable_;   % hash table: "variable name" -> [hash table: "list index" of handle list -> data]
  end % properties

  methods

    function [this] = ExtendedHashtable(arg)
      % Copy constructor
      if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
      %

      % use containers.Map object for hashtable "varname" -> Hashtable
      this.dataTable_   = Hashtable(); %containers.Map('KeyType','char','ValueType','any'); %Hashtable();
      %this.dataTable_ = containers.Map('KeyType', 'char', 'ValueType', 'any');  % seems to be more robust in very rare special cases
    end

    function Set(this, ename, evalue, ehandle)
    % SET
    %
    %  ename         - variable name (e.g. 'A')
    %  evalue        - value of name (e.g. MueMat Operator)
    %  ehandle       - handle of generating factory (e.g. SaPFactory object)
        if ~varexist('ehandle')
            ehandle = 0;
        end

        % create ename if it doesn't exist already
        if ~this.dataTable_.isKey(ename)
            this.dataTable_.Set(ename,Hashtable());
        end

        %this.dataTable_.Get(ename).Set(ehandle,evalue);
        dataTable = this.dataTable_.Get(ename);
        dataTable.Set(ehandle,evalue);
    end

    function [evalue] = Get(this, ename, ehandle)
    % GET
    %
    %  ename         - variable name (e.g. 'A')
    %  ehandle       - handle of generating factory (e.g. SaPFactory object)
    %  evalue        - return: value of name
        if ~varexist('ehandle')
            ehandle = 0;
        end

        % check if ename for ehandle is existent
        if ~this.dataTable_.isKey(ename)
           error('key %s not existent',ename);
        elseif ~this.dataTable_.Get(ename).isKey(ehandle)
            % should only happen if Needs cannot find a unique answer to a
            % call of "get".
            error('key %s not existent for specific handle',ename);
        end

        evalue = this.dataTable_.Get(ename).Get(ehandle);
    end

    function Remove(this, ename, ehandle)
    % REMOVE
    %
    %  ename         - variable name (e.g. 'A')
    %  evalue        - value of name
    %  ehandle       - handle of generating factory (e.g. SaPFactory object)
    %
    % removes variable 'ename', that is generated by ehandle factory
        if ~varexist('ehandle')
            ehandle = 0;
        end

        if ~this.isKey(ename, ehandle)
            error('%s not available',ename);
        end

        % 2) ename exists: ok
        %    check if ehandle already known for ename
        dataTable   = this.dataTable_.Get(ename);

        if dataTable.isKey(ehandle)
            dataTable.Remove(ehandle);
        end

        if length(dataTable) == 0
            % remove ename
            this.dataTable_.Remove(ename);
        end
    end


    function [ToF] = isKey(this, ename, ehandle)
    % isKey
    %
    %  ename         - variable name (e.g. 'A')
    %  evalue        - value of name
    %  ehandle       - handle of generating factory (e.g. SaPFactory object)
    %
        if ~varexist('ehandle')
             ehandle = 0;
        end

        % 1) check if ename is existent
        if ~this.dataTable_.isKey(ename)
            ToF = false;
            return;
        end

        % 2) ename exists: ok
        %    check if ehandle already known for ename
        if this.dataTable_.Get(ename).isKey(ehandle), ToF = true; return; end;

        % 3) ename not existent
        ToF = false;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    % helper functions for print, etc.

    function [keys] = keys(this)
    % keys
    %
    % returns cell list with all available keys (= varnames)
        keys = this.dataTable_.GetKeyList();
    end

    function [handles] = handles(this, ename)
    % handles
    %
    %  ename         - variable name (e.g. 'A')
    %
    % returns cell list with all available MATLAB handles for variable
    % 'ename'

        if ~this.dataTable_.isKey(ename)
            error('error');
        end

        dataTable = this.dataTable_.Get(ename);
        handles = dataTable.GetKeyList();
    end

  end % public methods

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

  end % protected methods

end %classdef