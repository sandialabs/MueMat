classdef Poisson_AMG < CopiableHandle
% This class provides a small facade class for a multigrid algorithm that can be used for elasticity problems

    methods (Static = true)
        function [MgHierarchy] = Setup(Amat,Sm,maxLevels)
            %Elasticity_AMG Setup
            %
            %   SYNTAX [MgHierarchy] = Setup(Amat,nsp,Sm,maxLevels)
            %
            %    Amat      - MueMat Operator with problem matrix
            %    Sm        - optional: Smoother (default: Smoother('GaussSeidel'))
            %    maxLevels - number of max. multigrid levels (default: 5)
            if ~varexist('Amat')
                error('Parameter Amat missing');
            end
	        if ~isa(Amat,'Operator')
		        error('Amat must be a MueMat Operator!');
	        end

            spacedim = Amat.GetRowMap().ConstBlkSize();

            if ~varexist('Sm')
                Sm = Smoother('GaussSeidel'); % default: 2 GaussSeidel iterations (no damping)
            end

            if ~varexist('maxLevels')
                maxLevels = 3;
            end

            Finest = Level;
            Finest.Set('A', Amat);

            Pfact = SaPFactory(TentativePFactory());	% standard smoothed aggregation AMG for symmetric problems
            Rfact = TransPFactory();

            MgHierarchy = Hierarchy();
            MgHierarchy.SetLevel(Finest,1);
            MgHierarchy.SetMaxCoarseSize(0.1*Amat.GetRowMap().NDOFs());
            MgHierarchy.FillHierarchy(Pfact,Rfact,RAPFactory(),1,maxLevels);
            MgHierarchy.SetSmoothers(SmootherFactory(Sm),1,maxLevels);

        end
    end

end
