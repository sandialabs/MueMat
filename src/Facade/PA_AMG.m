classdef PA_AMG < CopiableHandle
% This class provides a small facade class for a plain aggregation multigrid algorithm
% with identity matrix as default nullspace

    methods (Static = true)
        function [MgHierarchy] = Setup(Amat,Sm,maxlevels)
            if ~varexist('Amat')
               error('Parameter Amat missing');
            end
            if ~varexist('Sm')
               Sm = Smoother('GaussSeidel'); % default: 2 GaussSeidel iterations (no damping)
            end
            if ~varexist('maxlevels')
               maxlevels = 3;
            end
            if ~isa(Amat,'Operator')
               Amat = Operator(Amat);   % attention: number of null space vectors?
               warning('number of Dofs per node supposed to be 1');
            end

            Finest = Level();
            Finest.KeepAll(false);
            Finest.Set('A',Amat);

            MgHierarchy = Hierarchy();
            MgHierarchy.SetLevel(Finest,1);
            MgHierarchy.SetMaxCoarseSize(1);
            MgHierarchy.FillHierarchy(TentativePFactory(), TransPFactory(), RAPFactory(), 1, maxlevels);
            MgHierarchy.SetSmoothers(SmootherFactory(Sm));
        end
    end

end