classdef CrossFactory < CopiableHandle
% Handles parameter specifications that affect several factories
% Normally,
% only Hierarchy.m creates one though several factories use the static functions.
% Potential situations include
%
%    1) Needing a factory to save data (which is not normally saved)
%       so that a class, function, or factory can use it. Some examples are ...
%
%       a) aggregates saved for plotting.
%       b) tentative prolongators saved by prolongator factory for use within
%          restrictor factory.
%       c) prolongator sparsity pattern saved by prolongator factory so that a
%          related linear solve can reuse the pattern during its AMG setup.
%
%    2) Specify output or debugging level for all factories..
%
% One important point is that 'Needs' can come from users or be generated
% internally. For example, saving aggregates for plotting is a user-need while
% keeping aggregates for a restriction factory is an internal-need.  In the
% later case, the user may have no idea that a particular algorithm requires
% aggregates to be saved.
%
% See also CrossFactorySpecTest.m containing detailed documentation on its use.

  properties (Access = private)
     Specifications_   % structure holding needs (key/value) where
                       % Specifications_.key = value.
  end
  methods
     function [this] = CrossFactory(obj)
       % Copy constructor
       if nargin == 1 && isa(obj, class(this)), this.Copy_(obj,[]); return; end
       %

     % constructor
        this.Specifications_ = [];
     end
     function AddSpecs(this,Needs)
     % Adds new specifications to the list.
        this.Specifications_= CrossFactory.MergeNeeds(Needs,...
                                               this.Specifications_);
     end
     function [z] = GetSpecifications(this)
     % Gets entire specification list.
        z = this.Specifications_;
     end
     function [status] = TrueOrFalse(this,key,LevelId)
     % Checks if key or strcat('Persistent',key) is in specification list.
     % If in the list and the associated value is not of the form
     % a:b (where a & b are numbers), then it returns true. If the value has
     % the a:b form, then it returns true only if  a <= LevelId <= b.

       Persistentkey = strcat('Persistent',key);
       if isfield(this.Specifications_,Persistentkey), key=Persistentkey; end
       status = false;
       if isfield(this.Specifications_,key)
          if ~ischar(key), status = true;
          else
             eval(sprintf('str = this.Specifications_.%s;\n',key));
             if ischar(str)
                [range,count]  = sscanf(str,'%d:%d');
                if count ~= 2, status = true;
                elseif (LevelId>=range(1))&&(LevelId<=range(2)),status=true;end
             else
                status = false;
             end
          end
       end
     end
     function [Level] = MakeRequests(this, Level)
         if isempty(this.Specifications_), allfields = [];
         else allfields = fieldnames(this.Specifications_); end;
         id = Level.GetLevelId();
         for i=1:length(allfields)
           if this.TrueOrFalse(allfields{i},id)
              Level.Request(allfields{i})
              Level.Request(allfields{i})
              Level.Request(allfields{i})
              Level.Request(allfields{i})
              Level.Request(allfields{i})
              Level.Request(allfields{i})
           end
         end
     end
  end
  methods (Static = true)
     function [BNeeds] = MergeNeeds(ANeeds,BNeeds)
     % Merges two specification lists
     % A merge means something funny when both lists have the same key.
     % In particular,
     %     if value1, value2 are numeric | max(value1,value2)
     %     --------------------------------------------------
     %     if value1, value2 are         |   value1
     %        identical strings          |
     %     --------------------------------------------------
     %     if value1, value2 are strings | max(value1,value2)
     %        with 1 numerical value each|
     %     --------------------------------------------------
     %     if value1/value2 have a:b form| min(value1.a,value2.b):max(value1.a,value2.b):
     %        where a and b are numeric  |
     %     --------------------------------------------------
     %     if value1, value2 are         | strcat(value1,value2)
     %        non-numeric strings        |
     %     --------------------------------------------------
     %     else                          | value2
     %     --------------------------------------------------
        if ~isempty(ANeeds), list = fieldnames(ANeeds);
        else   list = []; end

        for i=1:length(list)
           if ~isfield(BNeeds,list(i)),
             BNeeds.(deblank(list{i}))=ANeeds.(deblank(list{i}));
           elseif isnumeric(ANeeds.(list{i})) && isnumeric(BNeeds.(list{i})),
                 BNeeds.(list{i})= max(ANeeds.(list{i}),BNeeds.(list{i}));
           elseif ischar(ANeeds.(list{i})) && ischar(BNeeds.(list{i})),
             [Alowhigh,Acount]=sscanf(ANeeds.(deblank(list{i})),'%d:%d');
             [Blowhigh,Bcount]=sscanf(BNeeds.(deblank(list{i})),'%d:%d');
             if (Acount == 2) && (Bcount == 2),
               BNeeds.(deblank(list{i})) = ...
                   sprintf('%d:%d',min(Blowhigh(1),Alowhigh(1)),...
                                   max(Blowhigh(2),Alowhigh(2)));
             else
                [Anum,Acount]=sscanf(ANeeds.(deblank(list{i})),'%d');
                [Bnum,Bcount]=sscanf(BNeeds.(deblank(list{i})),'%d');
                if (Acount == 1) && (Bcount == 1),
                  BNeeds.(deblank(list{i})) = sprintf('%d', max(Anum,Bnum));
                else
                 if ~strcmp(ANeeds.(deblank(list{i})),BNeeds.(deblank(list{i}))),
                  BNeeds.(deblank(list{i}))=strcat(ANeeds.(deblank(list{i})),...
                                                   BNeeds.(deblank(list{i})));
                 end
                end
             end % if (Acount == 2) && (Bcount == 2),
           end % if ~isfield(BNeeds,list(i)),
        end %for
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
