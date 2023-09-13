classdef TipSmoother < Smoother
  % TipSmoother class
  %
  % This smoother depends on Level.GetTips(). In particular,
  %
  %     Level.GetTips() is empty   :  Standard smoothing is done as 
  %                                   described in the Smoother class
  %     Level.GetTips() is nonempty:  A block variant of the standard 
  %                                   smoothing (as described in Smoother
  %                                   class) where the blocks are defined
  %                                   according to the Collection 
  %                                   associated with Level.GetTips().
  %
  
  properties (Access = private)
    
    TipDofs_

  end
  
  methods
    function [this] = TipSmoother(varargin)
      %SMOOTHER Constructor
      %
      %   SYNTAX   TipSmoother(string, nIts, Omega, diagonalView);
      %
      %     string         - smoother name (string)
      %     nIts           - number of iteration (integer, optional, default=2)
      %     Omega          - damping factor (double, optional, default=1)
      %     diagonalView   - Block or point smoother (optional)
      %
      %   EXAMPLE
      %
      %     Smoo = TipSmoother('GaussSeidel', 2, 1);
      
      % Copy Constructor code line unneeded here. It's done directly by this@Smoother();
           
      this@Smoother(varargin{:});

      if ~(nargin == 1 && isa(varargin{1}, class(this))) % if it is not the copy constructor, then:
        
        % This smoother use a special view. The view can be named by
        % the user using the input argument 'diagonalView'
        % By default, the name of the view is 'tips'
        if nargin < 4, this.SetDiagonalView('tips'); end; % Note: this.SetDiagonalView() will set isSetup flag to false but not a problem here.
      
      end

    end
    
    function SetNeeds(this, CurrentLevel)
      % Obtain any cross factory specifications
      %CurrentLevel.Request('Tips'); % call request from outside!!!
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Setup(this, Level)
      %SETUP Run the setup phase of the smoother.
      % See also: SmootherPrototype.Setup
      % 
      %   SYNTAX   Amat = obj.Setup(Level);
      % 
      %     Level - level of the MG hierachy (Level)
      
      if this.isSetup(), Level.Release('Tips'); return; end 

      tips = Level.Get('Tips');
      fprintf('# diagonal blocks in smoother = %d\n',tips.NSubsets);
      A = Level.Get('A');
      viewName = this.GetDiagonalView();
      
      if ~ismember(viewName, A.GetViewList())
         
         % Create a view 'tips' associated with the collection 'Level.GetTips()'
         A.CreateView(viewName, A.GetRowMap(), A.GetColMap(), A.GetApply());
         
         % Setup the special diagonal to store the collection
         previousView = A.SwitchToView(viewName);
         BlkDiag = A.GetDiagonal(tips);
         if isempty(BlkDiag.GetApplyInverse())
            FactorBlkDiag(BlkDiag);
         end
         A.SwitchToView(previousView);
         
      else
         warning(['TipSmoother: a view named ' this.GetDiagonalView() ' already exist.']);
         % TODO: check if the existing view is the same as the view I want to use.
      end

      % Real setup phase of the underlying smoother
      Setup@Smoother(this, Level);
      
      % Release
      Level.Release('Tips');

    end % function
    
  end % methods
  
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
  
end % class
