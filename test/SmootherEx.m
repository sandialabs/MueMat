classdef SmootherEx < Smoother 
    % Jacobi smoother with automatic choice of "optimal" damping parameter
    %                              2
    %     omega = --------------------------------------
    %             ( 1 + 1/coarseningrate) * lambda_{max}
    %
    % and lambda_{max} the maximum eigenvalue of the symmetric part of
    % D^{inv} A
    % The coarsening rate is calculated using the aggregation information
    % This smoother needs access to Level.Get('Aggregates')!
    % TODO: introduce "switch" to allow setting coarseing rate by hand
    % (e.g. for geometric multigrid).
    
    properties (Access = private)
        % Parameters of the smoother
        CoarseningRate_ = 3;       % coarsening rate
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [this] = SmootherEx(nIts, diagonalView)
            if ~varexist('diagonalView'), diagonalView = 'current'; end;
            
            this = this@Smoother('Jacobi', nIts, -1.0, diagonalView);
            
            % Copy Constructor
            if nargin == 1 && isa(nIts, class(this)), this.Copy_(nIts,[]); return; end
            %
            
            if varexist('nIts'),  this.SetNIts(nIts); end
            if varexist('diagonalView'), this.SetDiagonalView(diagonalView); end;
                       
            this.SetType('Jacobi (ext)');
        end% Smoother()
        
        function SetNeeds(this, Level)
            % Obtain any cross factory specifications
            Level.Request('Aggregates');
        end
   
        function [CoarseningRate] = GetCoarseningRate(this)
            %GETCOARSENINGRATE Get the coarsening rate
            %
            %   SYNTAX   CoarseningRate = obj.GetCoarseningRate();            
            CoarseningRate = this.CoarseningRate_;
        end
        
        function CopyParameters(this, src)
            %COPYPARAMETERS Copy the parameters of another smoother prototype.
            % See also: SmootherPrototype.CopyParameters
            %
            %   SYNTAX   obj.CopyParameters(src);
            %
            %     src - Object (SmootherPrototype)
            CopyParameters@Smoother(this,src);
            this.CoarseningRate_ = src.GetCoarseningRate();
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Setup(this, Level, Specs)
            %SETUP Run the setup phase of the smoother.
            % See also: SmootherPrototype.Setup, Smoother.SetupDiagonal, Smoother.Setup
            %
            %   SYNTAX   Amat = obj.Setup(Level, Specs);
            %
            %     Level - level of the MG hierachy (Level)
            %     Specs - specifications (CrossFactory)
            
            % This Smoother object might have an alternative definition
            % of the block diagonal or an alternative definition of 'block'
            % as well as a different damping parameter and number of iterations
            % associated with Smoother.Apply() depending on factory specifications.
      
            % TODO: check me!
            %             if this.isSetup(), 
            %                 Level.Release('Aggregates'); % release aggregates
            %                 return; 
            %             end 
            
            Setup@Smoother(this,Level);
            
            Amat = this.GetA();
            
            this.SetupDiagonal(); % 'view' cuisine
           
            % 1 ) calculate preconditioned matrix 
            BlkDiag = Amat.GetDiagonal([], this.GetDiagonalView());   % diagonalView should be "point"-based?
            DinvA = BlkDiag \ Amat;
            
            % 2) get symmetric part of D^{-1}A
            AA = DinvA.GetMatrixData();
            Asym = 0.5*(AA+AA');
            
            % 3) calculate biggest EW of Asym (approximatively)
            D = eig(full(Asym));
            lmax = max(D);
            
            % 4) get coarsening rates from aggregates
            AggInfo = Level.Get('Aggregates'); Level.Release('Aggregates');
            this.CoarseningRate_ = size(AggInfo.NodesInAgg,2)/size(AggInfo.NodesInAgg,1);
            
            % 5) approximate smallest eigenvalue for multigrid method
            %     
            %       \lambda_{min} \approx \frac{\lambda_{max}}{CoarseningRate}
            lmin = lmax/this.CoarseningRate_;
            
            % 6) calculate optimal omega
            %
            %     \omega_{opt} = \frac{2}{\lambda_{min} + \lambda_{max}}
            %            \approx
            %            \frac{2}{\frac{\lambda_{max}{CoarseningRate} + \lambda_{max}}
            %            = 2 \frac{1}{\bigl(\frac{1}{CoarseningRate}+1\bigr) \lambda_{max}}
            %
            %this.Omega_ = 2/(lmin+lmax); %3/(2*lmax);% 1/(lmin+lmax); % why not 2/(lmin+lmax???)  ->  
            this.SetOmega(2/(lmin+lmax));
            
            this.SetIsSetup(true);
            
            
        end
        
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

% TODO: cleanup comments
