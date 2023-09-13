%COPIABLEHANDLE Superclass to add copy capability to Matlab handle classes
%
% How to copy an object:
%
%   F = Foo(); cpyF = F.Copy() or cpyF = Foo(F)
%
% How to implement a copy constructor and a copy method using this class:
%
%   1) Instead of deriving from Matlab handle class, derive from this class.
%   2) Create a protected method Copy_() which copy every properties of
%      your class and call copy constructor of base classes.
%      You can create your own method or use this generic method:
%
%      methods (Access = protected)
%       function Copy_(this, src, mc)
%         [cmd, data, mc] = this.CopyCmd_(src,mc);
%         eval(cmd);
%       end
%      end % methods
%
%      This generic method creates deep copy of your object properties using
%      the utility method CopyCmd_ defined on this file.
%
%   3) Add at the begining of your constructor:
%    if nargin == 1 && isa(a, class(this)), this.Copy_(arg,[]); return; end
%      where arg is the first argument of your constructor (add it if any).
%
%  See also: Foo.m

% Misc:
%   - If developpers don't want a deep copy of every properties, they can
%     create their own implementation of Copy_(). The use of CopyCmd_() is
%     optional.
%
%   - Using such class is paintful. I contacted Matlab support and
%     the only way to copy classes is to implement a copy() method
%     in every class.
%
%   - According to MATLAB technical support, copiable handle will
%     be added to MATLAB in the future:
%
%     "The need for a "Copyable" type of superclass has already been
%     forwarded to the developers and they are looking into adding some
%     feature like that for a future release of MATLAB. Unfortunately,
%     I do not have an exact timeline as to when this feature will be
%     incorporated. However, the developers are actively looking into
%     this."
%
% Implementation notes:
%   - the only way to access private properties of objects is to convert
%     them to a struct.
%   - the only way to assign private properties is to do it from a method
%     of the class. That's why I need a Copy_() method on each class.
%   - the only way to create object is to use the class constructor.
%   - metaclass are handle objects describing MATLAB classes (properties,
%     superclasses...).
%
classdef CopiableHandle < handle

  methods (Access = protected, Static = true)

    function [cmd, data, mc] = CopyCmd_(src,mc)
      % Internal method to create generic copy constructor.
      %
      % This method return a string 'cmd' packing matlab code which:
      %  - copy every properties of your class
      %    (by value or with a call to 'this.Copy()' for handle classes).
      %  - call recursively Copy_() for base classes.
      %
      if ~isa(src,'struct')
        mc = eval(['?',class(src)]);
        origWarn = warning();
        warning('off', 'MATLAB:structOnObject');
        data = struct(src); % convert src to struct to access private properties
        warning(origWarn);
      else
        if ~varexist('mc'), error('CopyCmd_ needs ''mc'''); end
        data = src; % data is already a struct
      end

      cmd='';

      % Copy properties
      % note: the Copy function of utils/ is called for non-objects
      % properties. This function is equivalent to 'cpy = src'.
      for i=1:size(mc.Properties,1)
        cmd = [cmd ...
          'this.' mc.Properties{i}.Name ' = Copy(data.' mc.Properties{i}.Name '); '];
      end

      % Call recursively Copy_() for base classes.
      for i=1:length(mc.SuperClasses)
        cmd = [cmd ...
          'this.Copy_@' mc.SuperClasses{i}.Name '(data, mc.SuperClasses{' num2str(i) '}); '];
      end

    end

  end

  methods

    function cpy = Copy(this)
      % Create a copy of the calling object.
      % This method is equivalent to cpy = obj.Copy()

      % Call the copy constructor: copy = Foo(this)
      Constructor = class(this);
      cpy = eval([Constructor '(this)']);
    end %function Copy

  end %methods

  methods (Access = protected)
    function [cmd, data, mc] = Copy_(this,src,mc)
      % Internal Copy_ method for ClassBase objects.

      % nothing to do here.
    end
   end

end

% Some documentations:
%
% http://stackoverflow.com/questions/247430/matlab-copy-constructor
% * example of function copy
%   => does not work for 'private properties'
%
% http://stackoverflow.com/questions/2388409/how-can-i-tell-how-much-memory-a-handle-object-uses-in-matlab
% * example of how to convert class to struct to get private properties + set warning on/off
%
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/257925
% * Q/A: How to copy a handle object ?
%   - fns = properties(src);
%   - save/load on disk :(
%   - meta-class
%  => some very good answers but function copy needs full access to all of the properties and so should be added as a class method.
%
% http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_oop/br8b90p.html
% * Documentation about Meta-class
%
% http://www.advancedmcode.org/object-oriented-programming-in-matlab.html#95
% * Documentation about Meta-class
% * example of shallow copy using Meta-class
%
% http://www.mathworks.fr/matlabcentral/fileexchange/22965-clone-handle-object-using-matlab-oop
