classdef Elasticity_AMG < CopiableHandle
% This class provides a small facade class for a multigrid algorithm that can be used for elasticity problems

    methods (Static = true)
        function [MgHierarchy] = Setup(Amat,nsp,Sm,maxLevels)
            %Elasticity_AMG Setup
            %
            %   SYNTAX [MgHierarchy] = Setup(Amat,nsp,Sm,maxLevels)
            %
            %    Amat      - MueMat Operator with problem matrix
            %    nsp       - nullspace (rigid body modes)
            %    Sm        - optional: Smoother (default: Smoother('GaussSeidel'))
            %    maxLevels - number of max. multigrid levels (default: 5)
            if ~varexist('Amat')
               error('Parameter Amat missing');
            end

            if ~varexist('nsp')
                error('missing nullspace information!');
            end

            spacedim = Amat.GetRowMap().ConstBlkSize();
            if ((spacedim == 2 && size(nsp,2) ~=3) || ...
                (spacedim == 3 && size(nsp,2) ~=6))
                 error('space dimension and number of nullspace vectors for elasticity problems (rigid body modes) do not fit together');
            end

            if ~varexist('Sm')
                Sm = Smoother('GaussSeidel'); % default: 2 GaussSeidel iterations (no damping)
            end

            if ~varexist('maxLevels')
                maxLevels = 5;
            end
            if ~isa(Amat,'Operator')
               Amat = Operator(Amat);   % attention: number of null space vectors?
               warning('number of Dofs per node supposed to be 1');
            end

            Finest = Level;
            Finest.Set('A', Amat);
            Finest.Set('NullSpace', nsp);

            Pfact = SaPFactory();
            Rfact = TransPFactory();



            MgHierarchy = Hierarchy();
            MgHierarchy.SetLevel(Finest,1);
            MgHierarchy.SetMaxCoarseSize(0.1*Amat.GetRowMap().NDOFs());
            MgHierarchy.FillHierarchy(Pfact,Rfact,RAPFactory(),1,maxLevels);
            MgHierarchy.SetSmoothers(SmootherFactory(Sm),1,maxLevels);

        end
    end

end