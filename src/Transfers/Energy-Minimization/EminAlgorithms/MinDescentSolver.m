classdef MinDescentSolver < VerboseObject
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

        % StepLengthStrategy_ can be constant, adaptive, global
        StepLengthStrategy_  = [];

        % very special/experimental/complicated stuff
        % use the same number of minimization iterations for restrictor
        % than for prolongator (and same step lengths, too)
        % this is the only place i need to know if algorithm is currently
        % in prolongation or restriction mode (prolongation_mode parameter)
        % :-(
        % for global step length strategy, true seems to be reasonable
        % for adaptive step length strategy, false is more robust...
        SimpleRestriction_ = true;
    end

    methods
        function this = MinDescentSolver(nIts, StepLength, StepLengthStrategy)
            if nargin == 1 && isa(nIts, class(this)), this.Copy_(nIts,[]); return; end
            if varexist('nIts'), this.nIts_ = nIts; end;
            if varexist('StepLength'), this.StepLength_ = StepLength; end;

            if varexist('StepLengthStrategy')
                if strcmp(StepLengthStrategy, 'constant') || ...
                        strcmp(StepLengthStrategy, 'adaptive') || ...
                        strcmp(StepLengthStrategy, 'global')   || ...
                        strcmp(StepLengthStrategy, 'local')
                    this.StepLengthStrategy_ = StepLengthStrategy;
                else
                    error('unknown StepLengthStrategy %s, must be either constant, adaptive, global or local',StepLengthStrategy);
                end
            else
                this.StepLengthStrategy_ = 'global'; % default
                this.SimpleRestriction_ = true;
            end

            if strcmp(this.StepLengthStrategy_,'adaptive')
                this.SimpleRestriction_ = false;
            end
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

        function [ToF] = SetSimpleRestriction(this,ToF)
            if varexist('ToF'),
                ToFold = this.SimpleRestriction_;
                this.SimpleRestriction_ = ToF;
                ToF = ToFold;
            else ToF = this.SimpleRestriction_; end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
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

            switch lower(this.StepLengthStrategy_)
                case 'constant'
                    [PP] = this.MinDescentConstant(A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode);
                case 'adaptive'
                    [PP] = this.MinDescentAdaptive(A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode);
                case 'global'
                    [PP] = this.MinDescentGlobal(A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode);
                otherwise
                    error('unknown StepLengthStrategy: %s, must be either constant, adaptive or global',this.StepLengthStratgy_);
            end

            P = Operator(PP, P0.GetRowMap(), P0.GetColMap(), P0.GetApply());
        end
    end

    methods (Access = private)
        function [PP] = MinDescentConstant(this, A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode)
            AA = A.GetMatrixData();
            PP = P0.GetMatrixData();
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            Bzero = FineLevel.Get('NullSpace');
            FineLevel.Release('NullSpace');

            for k=1:this.nIts_
                dir = -ddd * AA * PP;
%
                if ~isempty(this.SatisfyConstraintsFunc_),
                    dir2 = this.SatisfyConstraintsFunc_(dir, [], CoarseLevel, Pattern, B);
                end

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
        end

        function [PP] = MinDescentAdaptive(this, A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode);

            % Input
            AA = A.GetMatrixData();
            PP = P0.GetMatrixData();
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            Bzero = FineLevel.Get('NullSpace');
            FineLevel.Release('NullSpace');

            % working variables
            steplength = this.StepLength_;
            lastAPnorm = inf;

            % for simple_restriction
            CoarseLevel.Keep('MinDescentAdaptive_MinIter',this);
            CoarseLevel.Keep('MinDescentAdaptive_StepLengthReduction',this);

            for k=1:this.nIts_
                dir = -ddd * AA * PP;
                if ~isempty(this.SatisfyConstraintsFunc_),
                    dir2 = this.SatisfyConstraintsFunc_(dir, [], CoarseLevel, Pattern, B);
                end

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

                %% step length control routine
                if prolongation_mode || ~this.SimpleRestriction_
                    % we are in "prolongation mode", just reduce step
                    % length size if needed
                    numStepLengthShortening = 0;
                    curAPnorm = norm(AA*(PP+steplength*dir),'fro');
                    while (curAPnorm > lastAPnorm)
                        steplength = 0.5*steplength;
                        numStepLengthShortening = numStepLengthShortening + 1;
                        curAPnorm = norm(AA*(PP+steplength*dir),'fro');
                        if(steplength < 1e-10)
                            break;
                        end;
                    end;
                    lastAPnorm = curAPnorm;
                    CoarseLevel.Set('MinDescentAdaptive_StepLengthReduction',numStepLengthShortening, this);
                else
                    % we are in "restriction mode": apply
                    % MinDescentAdaptive_MinIter iterations with full step
                    % length (simple restriction improvement)
                    if CoarseLevel.IsAvailable('MinDescentAdaptive_MinIter',this)
                        if  k>=CoarseLevel.Get('MinDescentAdaptive_MinIter',this)
                            % after MinIter full steps reduce step length
                            numStepLengthShortening = CoarseLevel.Get('MinDescentAdaptive_StepLengthReduction',this);
                            steplength = (0.5^numStepLengthShortening) * steplength;
                        end
                    end
                end

                PP = PP + steplength * dir;

                %% step length break condition & output
                if ~this.SimpleRestriction_
                    fprintf('MinDescent(adaptive) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,this.StepLength_,norm(AA*PP,'fro'));
                    if(steplength < 1e-10)
                        break;
                    end
                else

                    if prolongation_mode
                        CoarseLevel.Set('MinDescentAdaptive_MinIter', this.nIts_, this);
                        if (steplength < 1e-10)
                            CoarseLevel.Set('MinDescentAdaptive_MinIter', k, this);
                            break;
                        end
                        %if this.GetOutputLevel() > 5
                        fprintf('MinDescent(adaptive) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,this.StepLength_,norm(AA*PP,'fro'));
                        %end

                    else
                        if(steplength < 1e-10)
                            break;
                        end
                        %if this.GetOutputLevel() > 5
                        fprintf('MinDescent(adaptive) iter (restriction ): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,steplength,norm(AA*PP,'fro'));
                        %end
                    end
                end
            end
            if ~prolongation_mode
                % cleanup simple_restriction variables
                CoarseLevel.Delete('MinDescentAdaptive_MinIter',this);
                CoarseLevel.Delete('MinDescentAdaptive_StepLengthReduction',this);
            end
        end

        function [PP] = MinDescentGlobal(this, A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode);
            % Input
            AA = A.GetMatrixData();
            PP = P0.GetMatrixData();
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            Bzero = FineLevel.Get('NullSpace');
            FineLevel.Release('NullSpace');

            % working variables
            lastAPnorm = inf;

            % for simple_restriction
            CoarseLevel.Keep('MinDescentAdaptive_MinIter',this);
            CoarseLevel.Keep('MinDescentAdaptive_StepLengthReduction',this);

            for k=1:this.nIts_
                dir = -ddd * AA * PP;
                if ~isempty(this.SatisfyConstraintsFunc_),
                    dir2 = this.SatisfyConstraintsFunc_(dir, [], CoarseLevel, Pattern, B);
                end

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

                %% recalculate global steplength within every iteration
                % similar idea as in PG-AMG
                steplength = full(sum(diag(PP'*(AA')*AA*ddd*AA*PP)) / sum(diag(PP'*(AA')*ddd*(AA')*AA*ddd*AA*PP)));

                %% step length control routine
                if prolongation_mode || ~this.SimpleRestriction_
                    % we are in "prolongation mode", just reduce step
                    % length size if needed
                    numStepLengthShortening = 0;
                    curAPnorm = norm(AA*(PP+steplength*dir),'fro');
                    while (curAPnorm > lastAPnorm)
                        steplength = 0.5*steplength;
                        numStepLengthShortening = numStepLengthShortening + 1;
                        curAPnorm = norm(AA*(PP+steplength*dir),'fro');
                        if(steplength < 1e-10)
                            break;
                        end;
                    end;
                    lastAPnorm = curAPnorm;
                    CoarseLevel.Set('MinDescentAdaptive_StepLengthReduction',numStepLengthShortening, this);
                else
                    % we are in "restriction mode": apply
                    % MinDescentAdaptive_MinIter iterations with full step
                    % length (simple restriction improvement)
                    if CoarseLevel.IsAvailable('MinDescentAdaptive_MinIter',this)
                        if  k>=CoarseLevel.Get('MinDescentAdaptive_MinIter',this)
                            % after MinIter full steps reduce step length
                            numStepLengthShortening = CoarseLevel.Get('MinDescentAdaptive_StepLengthReduction',this);
                            steplength = (0.5^numStepLengthShortening) * steplength;
                        end
                    end
                end

                PP = PP + steplength * dir;

                %% step length break condition & output
                if ~this.SimpleRestriction_
                    fprintf('MinDescent(global) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,this.StepLength_,norm(AA*PP,'fro'));
                    if(steplength < 1e-10)
                        break;
                    end
                else

                    if prolongation_mode
                        CoarseLevel.Set('MinDescentAdaptive_MinIter', this.nIts_, this);
                        if (steplength < 1e-10)
                            CoarseLevel.Set('MinDescentAdaptive_MinIter', k, this);
                            break;
                        end
                        %if this.GetOutputLevel() > 5
                        fprintf('MinDescent(global) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,this.StepLength_,norm(AA*PP,'fro'));
                        %end

                    else
                        if(steplength < 1e-10)
                            break;
                        end
                        %if this.GetOutputLevel() > 5
                        fprintf('MinDescent(global) iter (restriction ): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',k,this.nIts_,steplength,norm(AA*PP,'fro'));
                        %end
                    end
                end
            end
            if ~prolongation_mode
                % cleanup simple_restriction variables
                CoarseLevel.Delete('MinDescentAdaptive_MinIter',this);
                CoarseLevel.Delete('MinDescentAdaptive_StepLengthReduction',this);
            end
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
            [cmd, dummy, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end

    end % methods
end