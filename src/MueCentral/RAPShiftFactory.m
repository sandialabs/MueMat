classdef RAPShiftFactory < RAPFactory

    % Special RAP Factory for indefinite Helmholtz systems
    % Factory which builds coarse discretizations via R*K*P-(alpha+i*beta)*(omega^2)*(R*M*P)
    %
    % Projects stiffness matrix (K) and mass matrix (M) onto coarse grid,
    % then recombines them to form A (we want to change the complex shift
    % at every level)
    % 
    % Assumes 'Afiltered' is the stiffness matrix K

    properties (Access = private)
        omega_  % frequency of the problem
        shifts_ % complex shifts at every level
    end

    methods
        function [this] = RAPShiftFactory()
            % constructor
	end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
        end

        function SetOmega(this, omega)
            this.omega_ = omega;
        end

        function SetShifts(this, shifts)
      	    this.shifts_ = shifts;
        end

        function flag = Build(this,FineLevel, CoarseLevel)
            flag = true;
            if ~CoarseLevel.IsAvailable('P')       fprintf('No Pmat??\n');             keyboard; end
            if ~CoarseLevel.IsAvailable('R')       fprintf('No Rmat??\n');             keyboard; end
	    if ~FineLevel.IsAvailable('A')         fprintf('No Amat??\n');             keyboard; end
            if ~FineLevel.IsAvailable('Afiltered') fprintf('No stiffness matrix??\n'); keyboard; end
            if ~FineLevel.IsAvailable('M')         fprintf('No mass matrix??\n');      keyboard; end
            Pmat   = CoarseLevel.Get('P');
            Rmat   = CoarseLevel.Get('R');
            omega  = this.omega_;
            shift  = this.shifts_(GetLevelId(CoarseLevel));
            Kmat   = FineLevel.Get('Afiltered');
            Mmat   = FineLevel.Get('M');
            Kc=Rmat*(Kmat*Pmat);
            Mc=Rmat*(Mmat*Pmat);
            CoarseLevel.Set('Afiltered', Kc);
            CoarseLevel.Set('M', Mc);
            CoarseLevel.Set('A', Kc-shift*(omega^2)*Mc); 
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
