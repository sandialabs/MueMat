classdef MultiVector < CopiableHandle
%MULTIVECTOR Class for supporting multivectors
% While not absolutely necessary in Matlab, this thin layer is
% needed in C++ to separate the multigrid algorithms from the underlying linear algebra library.
% This class should be used in non-compute intensive parts of MueMat.  In computational kernels, you should
% operate directly on the underlying data using the GetVectors() and SetVectors() methods.
% Vectors are stored column-major.

  properties (Access=private)
    data_            % vector data stored as column-major matrix
    numVectors_      % number of vectors
    vecLength_       % length of any single vector
  end %private properties

  methods

    %% Constructors

    function [this] = MultiVector(vecLength,numVectors,data,isSparse)
      %MULTIVECTOR Constructor.
      % Set number and length of vectors, and possibly load data or allocate memory.
      % Nonempty data takes precendence over specifying sizes, and there is no consistency check if you supply both.
      %
      %   SYNTAX   obj = MultiVector(vecLength, numVectors, data, isSparse);
      %
      %   EXAMPLE  obj = MultiVector()             % instance created, no memory allocated
      %            obj = MultiVector([],[],data)   % initialize from user-provided matrix (sparse or dense)
      %            obj = MultiVector(5,3)          % 3 dense vectors, each of length 5 and initialized to zero
      %            obj = MultiVector(5,3,true)     % 3 sparse vectors, each of length 5 (no memory allocated)
      %
      %     vecLength  - length of vectors (integer, optional, default=0)
      %     numVectors - number of vectors (integer, optional, default=0)
      %     data       - data of the MultiVector (MATLAB matrix, optional, default=[])
      %     isSparse   - sparse or dense (boolean, optional, default=false)

      % Copy Constructor
      if nargin == 1 && isa(vecLength, class(this)), this.Copy_(vecLength,[]); return; end

      if ~varexist('isSparse'), isSparse = false; end
      if nargin == 1
        error('MultiVector: specify both length and number of vectors, or data only');
      end
      if varexist('data') && ~isempty(data)
        this.data_ = data;
        this.numVectors_ = size(data,2);
        this.vecLength_ = size(data,1);
      elseif varexist('numVectors') && varexist('vecLength')
        if isSparse,  this.data_ = sparse(vecLength,numVectors);
        else          this.data_ = zeros(vecLength,numVectors);   end
        this.numVectors_ = numVectors;
        this.vecLength_ = vecLength;
      else
        this.data_ = [];
        this.numVectors_ = 0;
        this.vecLength_ = 0;
      end %if varexist(...
    end

    %% Mathematical operations

    function PutScalar(this,Scalar,ii)
      %PUTSCALAR Initialize all or part of multivector with the value 'Scalar'.
      % Note that if the entire multivector is initialized, the storage becomes dense.
      %
      %   SYNTAX   obj.PutScalar(Scalar, ii);
      %
      %     Scalar - value (double)
      %     ii     - iith vector of the MultiVector (integer, optional, default=1:numVectors)

      % FIXME check ii and Scalar
      if varexist('ii')
        this.data_(:,ii) = Scalar;
        %this.data_(:,ii) = 0;
      else
        this.data_ = ones(this.vecLength_,this.numVectors_) * Scalar;
        %[m,n] = size(this.data_);
        %if issparse(this.data_), this.data_ = sparse(m,n);
        %else                     this.data_ = zeros(m,n);  end
      end
    end

    function Random(this,ii)
      %RANDOM Randomize all or part of multivector.  This uses the default workspace stream.
      %
      %   SYNTAX   obj.Random(ii);
      %
      %     ii - iith vector of the MultiVector (integer, optional, default=1:numVectors)

      % FIXME check ii
      [m,n] = size(this.data_);
      if varexist('ii')
        if issparse(this.data_), this.data_(:,ii) = sprand(m,1,0.5);
        else                     this.data_(:,ii) = rand(m,1); end
      else
        if issparse(this.data_), this.data_ = sprand(m,n,0.5);
        else                     this.data_ = rand(m,n); end
      end
    end

    % methods we'll need
    %   add, subtract         ...            done
    %   norm                  ...            done
    %   transpose             ...            done
    %   multiply              ...            done
    %   component multiply    ...            done
    %   index
    %   is equal
    %   negative of           ...            done

    function [result] = plus(u,v)
      %PLUS Addition
      %
      %   SYNTAX   result = u + v;
      %
      %     u      - a MultiVector (or MATLAB matrix)
      %     v      - a MultiVector (or MATLAB matrix)
      %     result - a MultiVector
      [uleng,unumvec] = size(u);
      [vleng,vnumvec] = size(v);
      if unumvec ~= vnumvec || uleng ~= vleng
        error('MultiVector add: check the dimensions of your vectors')
      end
      if  isa(u,'MultiVector'), udata = u.GetVectors();
      else                      udata = u; end
      if  isa(v,'MultiVector'), vdata = v.GetVectors();
      else                      vdata = v; end
      result = MultiVector( [], [], udata+vdata);
    end

    function [result] = uminus(u)
      %UMINUS Unary minus
      %
      %   SYNTAX   result = -u;
      %
      %     u      - a MultiVector
      %     result - a MultiVector
      result = MultiVector([],[],-u.GetVectors());
    end

    function [result] = minus(u,v)
      %MINUS Subtraction.  Supports zero or one operand being a plain old Matlab vector.
      %
      %   SYNTAX   result = u - v;
      %
      %     u      - a MultiVector (or MATLAB matrix)
      %     v      - a MultiVector (or MATLAB matrix)
      %     result - a MultiVector
      [uleng,unumvec] = size(u);
      [vleng,vnumvec] = size(v);
      if unumvec ~= vnumvec || uleng ~= vleng
        error('MultiVector substract: check the dimensions of your vectors')
      end
      uIsMV = isa(u,'MultiVector');
      vIsMV = isa(v,'MultiVector');
      if uIsMV && vIsMV
        result = u + (-v);
      elseif uIsMV
        result = MultiVector( [], [], u.GetVectors()-v );
      else
        result = MultiVector( [], [], u-v.GetVectors() );
      end
    end

    function [result] = ctranspose(u)
      %CTRANSPOSE Transpose.
      %
      %   SYNTAX   result = u';
      %
      %     result - a MultiVector
      result = MultiVector([],[],u.GetVectors()');
    end

    function [result] = mtimes(u,v)
      %MTIMES Multiply.
      %
      %   SYNTAX   result = u * v;
      %
      %     u      - a MultiVector (or MATLAB matrix)
      %     v      - a MultiVector (or MATLAB matrix)
      %     result - a MultiVector
      [uleng,unumvec] = size(u);
      [vleng,vnumvec] = size(v);
      if unumvec ~= vleng
        error('MultiVector multiply: check the dimensions of your vectors')
      end
      if isa(u,'MultiVector'), udata = u.GetVectors();
      else                     udata = u;              end
      if isa(v,'MultiVector'), vdata = v.GetVectors();
      else                     vdata = v;              end
      result = zeros(vnumvec,1);
      for ii=1:vnumvec
        result(ii) = udata(ii,:) * vdata(:,ii);
      end
      % TODO should this return a MultiVector?
    end

    function [result] = times(u,v)
      %TIMES Component-wise multiply.
      %
      %   SYNTAX   result = u .* v;
      %
      %     u      - a MultiVector (or MATLAB matrix)
      %     v      - a MultiVector (or MATLAB matrix)
      %     result - a MultiVector
      [uleng,unumvec] = size(u);
      [vleng,vnumvec] = size(v);
      if unumvec ~= vnumvec || uleng ~= vleng
        error('MultiVector component multiply: check the dimensions of your vectors')
      end
      if isa(u,'MultiVector'), udata = u.GetVectors();
      else                     udata = u;              end
      if isa(v,'MultiVector'), vdata = v.GetVectors();
      else                     vdata = v;              end
      result = MultiVector([],[],udata .* vdata);
    end

    %% Attributes

    function [vecNorm] = norm(this,normType,ii)
      %NORM Return the norm of the iith vector.  Optionally accepts
      % any Matlab vector norm option (e.g., 1 or 'fro').
      %
      %   SYNTAX   vecNorm = obj.norm(normType, ii);
      %
      %     normType - type of norm (integer, optional, default=2)
      %     ii       - iith vector of the MultiVector (integer, optional, default=1:numVectors)
      %     vecNorm  - norm(s) of vector(s) (double or array of doubles)
      if ~varexist('normType'), normType = 2; end
      if varexist('ii'),
        this.CheckIndex(ii);
        vecNorm = norm(this.data_(:,ii),normType);
      else
        vecNorm = zeros(1,this.numVectors_);
        for ii=1:this.numVectors_,
          vecNorm(ii) = norm(this.data_(:,ii),normType);
        end
      end
    end

    function [numVectors] = NumVectors(this)
      %NUMVECTORS Return the number of vectors.
      %
      %   SYNTAX   numVectors = obj.NumVectors();
      %
      %     numVectors - number of vectors (integer)
      numVectors = this.numVectors_;
    end

    function [vecLength] = Length(this)
      %LENGTH Return the length of the MultiVector.
      %
      %   SYNTAX   vecLength = obj.Length();
      %
      %     vecLength  - length of vectors (integer)
      vecLength = this.vecLength_;
    end

    function [vecLength,numVectors] = size(u)
      %SIZE Returns length and number of vectors. Overloads the built-in size function.
      %
      %   SYNTAX   [vecLength, numVectors] = size(u);
      %
      %     u          - a MultiVector
      %     vecLength  - length of vectors (integer)
      %     numVectors - number of vectors (integer)
      vecLength = u.vecLength_;
      numVectors = u.numVectors_;
    end

    %% Access

    function SetVectors(this,data)
      %SETVECTORS Set MultiVector data.
      %
      %   SYNTAX   obj.SetVectors(data);
      %
      %     data - data of the MultiVector (MATLAB matrix)
      this.data_ = data;
      this.numVectors_ = size(data,2);
      this.vecLength_ = size(data,1);
    end

    function [data] = GetVectors(this,ii)
      %GETVECTORS Return MultiVector data.
      %
      %   SYNTAX   data = obj.GetVectors(ii);
      %
      %     ii   - iith vector of the MultiVector (integer, optional, default=1:numVectors)
      %     data - data of the MultiVector (MATLAB matrix)
      if varexist('ii')
        this.CheckIndex(ii);
        data = this.data_(:,ii);
      else
        data = this.data_;
      end
    end

    %function disp(this)
    %  % Print multivector to screen.  For a MultiVector "mv", this overloads
    %  % the action
    %  %   >> mv
    %  fprintf('    multivector\n');
    %  disp(this.data_)
    %end

  end %public methods

  methods (Access=private)

    function CheckIndex(this,ii)
      %CHECKINDEX
      %
      %   SYNTAX   obj.CheckIndex(ii);
      %
      %     ii - iith vector of the MultiVector (integer)
      if ii<1 || ii > this.numVectors_
        error( ['MultiVector.GetVectors: index ' num2str(ii) ...
               ' outside of range [1,' num2str(this.numVectors_) '].'] );
      end
    end

    function [udata] = GetData(this,u)
      %GETDATA
      %
      %   SYNTAX   udata = obj.GetData(u);
      %
      %     u     - a MultiVector
      %     udata - data of the MultiVector (MATLAB matrix)
      if isa(u,'MultiVector'), udata = u.GetVectors();
      else                     udata = u;              end
    end

  end %private methods

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

  end % method

end %class MultiVector
