classdef SteepestDescentSolver < VerboseObject
    properties (Access = private)
        nIts_ = 4;
        StepLength_ = 1.0;
        SatisfyConstraintsFunc_ = @EminPFactory.SatisfyConstraints

        % can be true or false
        % only optimization for non dirichlet Dofs?
        % true: we do not prescribe a wrong nullspace approximation for the
        % dirichlet bcs (default)
        % false: just "optimize" prolongator/restrictor for all Dofs even
        % though our optimal prolongation operators are disturbing the
        % solver at dirichlet bc dofs.
        IgnoreDirichletDofs_ = true;
    end

    methods
        function this = SteepestDescentSolver(nIts, StepLength)
            if nargin == 1 && isa(nIts, class(this)), this.Copy_(nIts,[]); return; end
            if varexist('nIts'), this.nIts_ = nIts; end;
            if varexist('StepLength'), this.StepLength_ = StepLength; end;
        end

        function SetNumIterations(this, nIts)
            this.nIts_ = nIts;
        end
        function SetStepLength(this, StepLength)
            this.StepLength_ = StepLength;
        end

        function SetNeeds(this,FineLevel, CoarseLevel)
            FineLevel.Request('NullSpace');
        end

        function [ToF] = IgnoreDirichletDofs(this, ToF)
            % ignore DiricheletDofs within descent minimization routine?
            % (get/set function)
            %
            %   SYNTAX ToF = IgnoreDirichletDofs(ToF)
            %
            %    ToF  - ignore DirichletDofs within descent minimization
            %           routine (true/false)

            if varexist('ToF') this.IgnoreDirichletDofs_ = ToF; end;
            ToF = this.IgnoreDirichletDofs_;
        end
        function P = Iterate(this, A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode)
            % P = Iterate(A, PatternFact, B, P0, CoarseLevel)
            %
            %    A            - MueMat operator (energy matrix)
            %    PatternFact  - pattern factory (for constraints)
            %    B            - constraint matrix
            %    P0           - initial prolongator (MueMat operator)
            %    CoarseLevel  - coarse level
            %    prolongation_mode - true, if EminPfactory is in
            %                        prolongation mode
            AA = A.GetMatrixData();
            PP = P0.GetMatrixData();
            %ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            Bzero = FineLevel.Get('NullSpace');
            FineLevel.Release('NullSpace');

            for k=1:this.nIts_
                dir = -2 * AA' * AA * PP;
                dir2 = this.SatisfyConstraintsFunc_(dir, [], CoarseLevel, Pattern);

                %% correct search direction (prolongator)
                % ignore dirichlet dofs, i.e. no disturbance for dirichlet
                % dofs?
                if this.IgnoreDirichletDofs_ == true

                    correctnsp = find(abs(AA*Bzero)<0.2);
                    dir(correctnsp,:) = dir2(correctnsp,:);
                    %dir(correctnsp) = dir2(correctnsp); % TODO error
                else
                    dir = dir2;
                end

                PP = PP + this.StepLength_ * dir;
                if this.GetOutputLevel() > 5
                    fprintf('MinDescent(constant) iter: %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.MinIter_,this.StepLength_,norm(AA*PP,'fro'));
                end
            end

            P = Operator(PP, P0.GetRowMap(), P0.GetColMap(), P0.GetApply());
        end
    end
end