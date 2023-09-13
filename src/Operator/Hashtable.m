classdef Hashtable < CopiableHandle
  %HASHTABLE A simple hashtable class
  %
  % Ex: t = Hashtable();
  %     t.Set('key1',magic(3));
  %     t.Set('key1',magic(4));
  %     t.GetKeyList()
  %     magic3 = t.Get('key1')
  %
  % This is actually not a hashtable at all. The implementation uses
  % an array (table_) and Get/Set methods look through the array to
  % find the matching key.
  %
  % Hashtable are only use for a very few number of keys, so we don't
  % care about performance here.

  properties (Access = private)
    table_ = struct('key_',[],'data_',[]);     % array of structures. Each structure has fields for key and data.
    size_ = 0 % number of entries in the array
  end % properties

  methods

    function [this] = Hashtable(arg)
      % Copy constructor
      if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
      %
    end

    function Set(this,key,data)
      % If an entry doesn't already exist for 'key', create a new one.
      ii = this.FindIndex(key);
      if ii == 0
        this.size_ = this.size_+1;
        this.table_(this.size_).key_  = key;
        this.table_(this.size_).data_ = data;
      else
        this.table_(ii).data_ = data;   % overwrite data
      end
    end

    function Remove(this,key)
       ii = this.FindIndex(key);
       if ii > 0
          this.table_(ii).key_ = -1;
          this.table_(ii).data_ = [];
          this.table_(ii) = [];
          this.size_ = this.size_-1;
       else
          fprintf('Hashtable::Remove: %s not found\n',key); keyboard;
       end
    end

    function data = Get(this,key)
      % Return [] is the key doesn't exist in the hashtable.
      ii = this.FindIndex(key);
      if ii == 0,
          data = [];
          fprintf('Hashtable::Get: %s not found\n',key); keyboard;
          return;
      end
      data = this.table_(ii).data_;
    end

    function [list] = GetKeyList(this)
      % Return a list (cell array) of the keys in the hashtable.
      if this.size_ == 0, list = []; return; end;
      for ii=1:this.size_
        list{ii} = this.table_(ii).key_;
      end
    end

    function [data] = GetDataList(this)
      % Return a list (cell array) of values in the hashtable.
      if this.size_ == 0, data = []; return; end;
      for ii=1:this.size_
        data{ii} = this.table_(ii).data_;
      end
    end

    function [ToF] = isKey(this,key)
      ToF = this.FindIndex(key) > 0;
    end

    function Print(this)
      for ii=1:this.size_
        this.table_(ii)
      end
    end

  end % public methods

  methods (Access = private)

    function [index] = FindIndex(this,key)
      % Return integer index into table corresponding to key "key".
      index=0;
      for ii=1:this.size_
        if isa(this.table_(ii).key_,'char')
            if strcmp(this.table_(ii).key_,key), index = ii; break; end
        else
            if this.table_(ii).key_ == key, index = ii; break; end
        end
      end
    end %FindIndex

  end % private methods

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
