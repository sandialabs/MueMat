
% Elasticity Test
function [ITERS, RESID, OC] = PoissonTest(dim, NLIST, stretch, EXPLIST, multimaterial, agg_rule)

if ~varexist('dim'),     dim = 2; fprintf('DIM = default = %d\n', dim); end

if ~varexist('NLIST')
  if dim == 1, NLIST = 10; end
  if dim == 2, NLIST = 50; end
    % if dim == 2, NLIST = 40:40:320; end
    if dim == 3, NLIST = 10;    end
    % if dim == 3, NLIST = 5:5:40;    end
end

if ~varexist('stretch'),       stretch = 1;       end % default = no stretching
if ~varexist('EXPLIST'),       EXPLIST = {'sa-nr', 'sa', 'emin', 'ml-sa', 'ml-sa-nr'}; end % default = run all
if ~varexist('multimaterial'), multimaterial = 0; end % default = steel only
if ~varexist('agg_rule'),      agg_rule = 0;      end % default = distLap

%numDesiredLevels = 10;
numDesiredLevels = 3;

%% Init environment
srand;
SetHomeDir
mue_include

%% Set test case parameters
% Parameters for 304 Stainless Steel
STEEL_MODULUS = 193e9;
POISSONS_RATIO = .305;

RUBBER_MODULUS = 1e8;

global numNS_EMIN;
numNS       = 1;
numNS_NOROT = 1;
numNS_EMIN  = 1;
wgts=1;


%% Set solver parameters
eminSteps = 2;
if stretch == 1
  AggDropTol = 1e-100;
else
  AggDropTol = 0.05;
end

%% Init output
ITERS = zeros(length(NLIST), 5);
OC = zeros(length(NLIST), 5);
NLEVELS = zeros(length(NLIST), 5);
LABEL = cell(1, 5);

for I = 1:length(NLIST)
    n = NLIST(I);
    global points_per_dim;
    points_per_dim = n;
    nstr = num2str(n);
    sstr = num2str(stretch);
    dimstr = num2str(dim);

    % Constant Elastic modulus
    if(multimaterial)
        % Steel and Rubber, approximately 50-50 split in the x direction
        %ELASTIC_MODULUS = @(x)((x(1)<(n+1)*stretch/2)*(STEEL_MODULUS-RUBBER_MODULUS)+RUBBER_MODULUS);

        % Steel and Rubber, sandwiched in ~6 layers in the y direction.
        %maxy = n+1;
        %ELASTIC_MODULUS = @(x)( (x(2)/maxy<.1 || (x(2)/maxy>.3&&x(2)/maxy<.5) || (x(2)/maxy>.7&&x(2)/maxy<.9))...
        %                      *(STEEL_MODULUS-RUBBER_MODULUS)+RUBBER_MODULUS);

        % Steel and Rubber, sandwiched in 3 layers in the y direction.
        maxy = n+1;
        miny = 2;
        ELASTIC_MODULUS = @(x)( ((x(2)-miny)/(maxy-miny)<.33333 || ((x(2)-miny)/(maxy-miny)>.66666))...
            *(STEEL_MODULUS-RUBBER_MODULUS)+RUBBER_MODULUS);

    else
        % 100% Steel
        ELASTIC_MODULUS = STEEL_MODULUS;
    end

    if dim == 2, sizestr = [nstr 'x' nstr]; else sizestr = [nstr 'x' nstr 'x' nstr]; end

    fprintf('***Running %s %d/%d***\n', sizestr, I, length(NLIST));

    % Build matrix, block it and grab parameters
    filename = [MUEMAT_ROOT_DIR '/data/Elas' dimstr 'D-' sizestr];
    if(stretch ~= 1) filename = [filename '.s.' sstr ];end
    if(multimaterial == 1) filename = [filename '.multmat' ];end
    filename = [filename '.mat'];

    fid = fopen(filename, 'r');
    % fid = -1; uncomment this line to always regenerate matrices
    fprintf('Generating matrix (%s)', filename);
    if dim == 1,
      Amat=BuildLaplace1D(n);
    elseif dim == 2
      Amat=BuildLaplace2D(n);
    end
    nullspace=ones(size(Amat.GetMatrixData,1),1);
    coords=[];


    fprintf('... done\n');

    %    Amat = Operator(Amat, 1, 1);
    RowMap = Amat.GetRowMap();
    ndofs = RowMap.NDOFs();

    % Create (normalized) RHS, run parameters
    rhs = rand(RowMap.NDOFs(), 1);
    rhs = rhs / norm(rhs);
    maxit = 50;
    tol = 1e-10;
    zeroGuess = zeros(RowMap.NDOFs(), 1);

    fprintf('CMS: Original Matrix NNZ = %d\n', nnz(Amat.GetMatrixData));

    global print_drop;

    if (ismember('sa-nr', EXPLIST))
        % Run SA w/o rotational mode
        J = 1;
        LABEL{J} = 'sa-norotate';
        [MgHierarchy, OC(I, J), NLEVELS(I, J)] = build_amg_hierarchy(dim, Amat.Copy(), nullspace(:, 1:numNS_NOROT), 'sa', numNS_NOROT, wgts, eminSteps, coords, AggDropTol, ELASTIC_MODULUS, agg_rule, numDesiredLevels);
        [sol, flag, relres, ITERS(I, J), RESID{I, J}] = pcg(Amat.GetMatrixData(), rhs, tol, maxit, ...
            @(b)MgHierarchy.Iterate(b, 1, zeroGuess, ALLZEROS, VCYCLE));
    end % sa-nr

    if (ismember('sa', EXPLIST))
        % Run Smoothed aggregation
        J = 2;
        print_drop = 1;
        LABEL{J} = 'sa';
        [MgHierarchy, OC(I, J), NLEVELS(I, J)] = build_amg_hierarchy(dim, Amat.Copy(), nullspace(:, 1:numNS), 'sa', numNS, wgts, eminSteps, coords, AggDropTol, ELASTIC_MODULUS, agg_rule, numDesiredLevels);
        [sol, flag, relres, ITERS(I, J), RESID{I, J}] = pcg(Amat.GetMatrixData(), rhs, tol, maxit, ...
            @(b)MgHierarchy.Iterate(b, 1, zeroGuess, ALLZEROS, VCYCLE));
        print_drop = 0;
    end % sa

    if (ismember('emin', EXPLIST))
        % Run Ray's special sauce
        J = 3;
        %LABEL{J} = ['svd-' num2str(numNS_EMIN)];
        LABEL{J} = 'emin';
        [MgHierarchy, OC(I, J), NLEVELS(I, J)] = build_amg_hierarchy(dim, Amat.Copy(), nullspace(:, 1:numNS), 'emin', numNS, wgts, eminSteps, coords, AggDropTol, ELASTIC_MODULUS, agg_rule, numDesiredLevels);
        [sol, flag, relres, ITERS(I, J), RESID{I, J}] = pcg(Amat.GetMatrixData(), rhs, tol, maxit, ...
            @(b)MgHierarchy.Iterate(b, 1, zeroGuess, ALLZEROS, VCYCLE));
    end % emin

    % Used by both ml-sa and ml-sa-nr:
    % Setting ML max levels
    % Make sure that the max number of levels used by ML is the same as MueLu
    nlindx = NLEVELS(I,:)~=0; nl = NLEVELS(I,nlindx);
    if (isempty(nl))
        fprintf('No MueMat run. Default numDesiredLevels for ML = %d\n', numDesiredLevels);
        ML_LEVELS = numDesiredLevels; % == default of MueMat
        ML_MAXCOARSE = 50; % == default of MueMat
        NLEVEL(I) = -1; % for output table -1 bcse we do not know how many levels are actually used by ML
    else
        if(~all(nl==nl(1))), fprintf('Error\n'); dbstack; keyboard; end % make sure all MueMat hierarchies have the same num of levels
        fprintf('ML numDesiredLevels == MueMat == %d\n', nl(1));
        ML_LEVELS = nl(1);     % force the number of Level
        ML_MAXCOARSE = 1;      % do not use this criteria
        NLEVEL(I) = nl(1);     % for output table
    end
    %

    % Run good ol' ML
    if 0%(ismember('ml-sa', EXPLIST))
        if ~exist('ml', 'file'), fprintf('This experiment requires MLMEX. Type ''doc MLInterface'' for more information.\n');
        else
            J = 4;
            LABEL{J} = 'ml-sa';
            if dim == 2
                ns = full([nullspace(:, 1);nullspace(:, 2);nullspace(:, 3)]);
                coords3 = [];
            else
                ns = full([nullspace(:, 1);nullspace(:, 2);nullspace(:, 3);nullspace(:, 4);nullspace(:, 5);nullspace(:, 6)]);
                coords3 = coords(:, 3);
            end

            [h, OC(I, J)] = ml('setup', Amat.GetMatrixData(), 'PDE equations', dim, ...
                'ML output', 10, 'smoother: type', 'symmetric Gauss-Seidel', ...
                'smoother: sweeps', 2, 'coarse: max size', ML_MAXCOARSE, 'max levels', ML_LEVELS, ...
                'coarse: type', 'Amesos-KLU', ...
                'aggregation: aux: enable', true, ...
                'aggregation: aux: threshold', AggDropTol, ...
                'x-coordinates', coords(:, 1), ...
                'y-coordinates', coords(:, 2), ...
                'z-coordinates', coords3, ...
                'null space: type', 'pre-computed', ...
                'null space: dimension', numNS, 'null space: vectors', ns);

            [sol, flag, relres, ITERS(I, J), RESID{I, J}] = pcg(Amat.GetMatrixData(), rhs, tol, maxit, ...
                @(b)ml(h, Amat.GetMatrixData(), b, 'krylov: type', 'fixed point', ...
                'krylov: output level', 0, ...
                'krylov: max iterations', 1, 'krylov: tolerance', 1e-100));

            %'aggregation: damping factor', 1e-100, ...% CMS
            %                 'aggregation: threshold', AggDropTol, ...

            ml('cleanup');
        end %if ~exist('ml', 'file')
    end % ml-sa

    if 0%(ismember('ml-sa-nr', EXPLIST))
        if ~exist('ml', 'file'), fprintf('This experiment requires MLMEX. Type ''doc MLInterface'' for more information.\n');
        else
            J = 5;
            LABEL{J} = 'ml-sa-nr';
            if dim == 2
                ns = full([nullspace(:, 1);nullspace(:, 2)]);
                coords3 = [];
            else
                ns = full([nullspace(:, 1);nullspace(:, 2);nullspace(:, 3)]);
                coords3 = coords(:, 3);
            end

            [h, OC(I, J)] = ml('setup', Amat.GetMatrixData(), 'PDE equations', dim, ...
                'ML output', 10, 'smoother: type', 'symmetric Gauss-Seidel', ...
                'smoother: sweeps', 2, 'coarse: max size', ML_MAXCOARSE, 'max levels', ML_LEVELS, ...
                'coarse: type', 'Amesos-KLU', ...
                'aggregation: aux: enable', true, ...
                'aggregation: aux: threshold', AggDropTol, ...
                'x-coordinates', coords(:, 1), ...
                'y-coordinates', coords(:, 2), ...
                'z-coordinates', coords3, ...
                'null space: type', 'pre-computed', ...
                'null space: dimension', numNS_NOROT, 'null space: vectors', ns);

            [sol, flag, relres, ITERS(I, J), RESID{I, J}] = pcg(Amat.GetMatrixData(), rhs, tol, maxit, ...
                @(b)ml(h, Amat.GetMatrixData(), b, 'krylov: type', 'fixed point', ...
                'krylov: output level', 0, ...
                'krylov: max iterations', 1, 'krylov: tolerance', 1e-100));

            % 'aggregation: damping factor', 1e-100, ...% CMS
            %                 'aggregation: threshold', AggDropTol, ...

            ml('cleanup');
        end %if ~exist('ml', 'file')
    end % ml-sa-nr

    % Output diagnostics
    fprintf('\nstatistics: ');
    for J = 1:length(LABEL), fprintf('%s %2d(%3.2f)   ', LABEL{J}, ITERS(I, J), OC(I, J)); end
    fprintf('\n');

    % Save results
    if(stretch == 1), save(sprintf('etest%dd-filter3_304.%s.mat', dim, cdate), 'LABEL', 'NLEVELS', 'NLIST', 'ITERS', 'RESID', 'OC');
    else save(sprintf('etest%dd-filter3_304.s%d.%s.mat', dim, stretch, cdate), 'stretch', 'LABEL', 'NLEVELS', 'NLIST', 'ITERS', 'RESID', 'OC'); end

    % vizAggs(MgHierarchy, 1, 8);

end % NLIST loop

fprintf('DIM = %d\n', dim);
fprintf('STRETCH = %f\n', stretch(1));

% Useful for printing
% TODO: code factorization with ElasticityTxtTable.m
fprintf('\n');
fprintf('        |  SA-NR  |ML-SA-NR ||    SA   |  ML-SA  ||   EMIN  |\n');
fprintf(' SZ LVL | ITS  OC | ITS  OC || ITS  OC | ITS  OC || ITS  OC |\n');
fprintf('-------------------------------------------------------------\n');
fprintf('%3d %2d  | %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f |\n',[NLIST',NLEVEL(:),...
    ITERS(:,1),OC(:,1),... % SA-NR
    ITERS(:,5),OC(:,5),... % ML-SA-NR
    ITERS(:,2),OC(:,2),... % SA
    ITERS(:,4),OC(:,4),... % ML-SA
    ITERS(:,3),OC(:,3)]'); % EMIN

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MgHierarchy, oc, NLevels] = build_amg_hierarchy(dim, Amat, nullspace, method_type, numNS, wgts, eminSteps, coords, AggDropTol, ELASTIC_MODULUS, agg_rule, numDesiredLevels)
%  Construct and populate finest level
Finest = Level();
Finest.Set('A', Amat);
%Finest.Set('NullSpace', nullspace);
%Finest.Set('xcoords', coords(:, 1));
%Finest.Set('ycoords', coords(:, 2));
%if dim == 3
%    Finest.Set('zcoords', coords(:, 3));
%end




% Set options
%aggoptions.AggTargetSize = 3;
%SAoptions.NoQR = 1;
%SAoptions.CNullMode = 'inject partialFN';
SAoptions.CNullMode = 'inject CMS';
%SAoptions.CNullMode = 'inject SVD';

% Build agglomeration factory
AmalgamateDropFact = CoalesceDropFactory();
%AmalgamateDropFact = CoalesceDoubleDropFactory();
%AmalgamateDropFact.SetDoublePostDropSpecifications(@SAdrop, AggDropTol);



% Weight function for the distance Laplacian to discourage
% aggregation across material boundaries

%if(agg_rule == 0),
%    fprintf('ElasticityTest: Using Distance Laplacian\n');
%    DISTANCE = @(x_, y_)(norm(x_-y_));
%elseif(agg_rule == 1),
%    fprintf('ElasticityTest: Using Material Laplacian\n');
%    DISTANCE = @(x_, y_)(sqrt((ELASTIC_MODULUS(x_)+ELASTIC_MODULUS(y_))/2));
%elseif(agg_rule == 2),
%    fprintf('ElasticityTest: Using Weighted Distance Laplacian\n');
%    BASE_EM = ELASTIC_MODULUS([0, 0, 0])^(1/3);
%    %DISTANCE = @(x_, y_)( norm(x_-y_) * ((ELASTIC_MODULUS(x_)*ELASTIC_MODULUS(y_))^(1/3)/(BASE_EM*BASE_EM)));
%    DISTANCE = @(x_, y_)( norm(x_-y_) *  sqrt((ELASTIC_MODULUS(x_)+ELASTIC_MODULUS(y_))/(2*BASE_EM)));
%else
%    fprintf('ElasticityTest: Unknown AggRule\n');
%    keyboard;
%end



% Amalgamate Amat to build auxillary coordinate Laplacian for aggregation
% Note: We do this *before* we do the dropping since we want the
% full distance Laplacian.  Also set the AuxMatrixFunc to
% regenerate (rather than coarsen) the distance Laplacian.
%AmalgMat = AmalgamateDropFact.Build_(Amat);
%AuxMat = BuildDistanceLaplacian(AmalgMat, coords);
%Finest.Set('AuxMatrix', AuxMat);
%Finest.Set('AuxMatrixFunc', @(A_, x_)BuildDistanceLaplacian(A_, x_, DISTANCE));
%AmalgamateDropFact.SetAName('AuxMatrix');

%fprintf('CMS: size(A) = %d size(AuxMatrix) = %d\n', size(Amat.GetMatrixData(), 1), size(AuxMat.GetMatrixData(), 1));


% Set dropping
AmalgamateDropFact.SetPreDropSpecifications(@SAdrop, AggDropTol);
%AmalgamateDropFact.SetPostDropSpecifications(@SAdrop, AggDropTol);

% Aggregation
% Note: Dropping should be done in ADF, not ML aggregation
AggFact = AggregationFactory();
AggFact.SetAlgorithm('ml'); % use the same aggregation as ML (need MLMEX)
AggFact.SetMLOptions('ML output', 10, 'aggregation: threshold', AggDropTol);

%global points_per_dim;
%AggFact.SetAlgorithm('order');
%AggFact.SetTargetSize({1, 3, 3});
%AggFact.SetTargetSize({1, 1, 3});
%AggFact.SetDomainPtsPerDim(points_per_dim*ones(dim, 1));


if(strcmp(method_type, 'sa')),
    % Standard Smoothed Aggregation

    %    filteredAuxMat = AmalgamateDropFact.Build_(AuxMat);
    %    filteredA = BuildDecoalescedFilteredMatrix(Amat, filteredAuxMat);
    %    Finest.Set('AfilteredDecoalesced', filteredA);

    Ptentfact = TentativePFactory(AmalgamateDropFact, AggFact);
    %Pfact = Ptentfact;

    Pfact     = SaPFactory(Ptentfact);
    %    Pfact.SetAForSmoothing('AfilteredDecoalesced');

    % This optim is not needed anymore:
    % % On the finest level, the block diagonal seems to have only zero
    % % outside of the main diagonal.
    % % Use the 'point' diagonal for faster computation.
    % % The same is done for the smoother.
    % PfactLvl1 = SaPFactory(Ptentfact,'point');
    % PfactLvl1.SetAForSmoothing('AfilteredDecoalesced');

elseif(strcmp(method_type, 'emin')),
    % Energy Minimization
    global numNS_EMIN
    SAoptions.NCoarseDofPerNode = numNS_EMIN;
    wgts = wgts(1:numNS);
    %fprintf('weights = [%s]\n', num2str(wgts));
    SAoptions.CoarseNullspaceWeightingMultiplier = wgts;
    SAoptions.DropConstraintsMethod = 'nullspace';
    SAoptions.FilteredAP = true;
    %SAoptions.PtentRootModifications = '4b'; % NOT SUPPORTED
    %SAoptions.SmoothNullspace = true;%CMS

    CNSFact = CoarseNSFactory();


    Pinitfact      = TentativePFactoryEx(AmalgamateDropFact, AggFact, CNSFact, SAoptions);
    PatFact       = AP_PatternFactory([], Pinitfact, SAoptions); % no filter, only options
    %PatFact        = DecoalescedGraph_PatternFactory([], AmalgamateDropFact, Pinitfact, SAoptions);

    FCSplittingFact = FCSplittingFactory(AggFact);
    ConstraintFact  = ConstraintFactory(PatFact, FCSplittingFact, SAoptions);


    Pfact          = EminPFactory(PatFact, ConstraintFact, CGEminSolver(eminSteps), Pinitfact, SAoptions);


    %% This optim is not needed anymore:
    %%  PfactLvl1      = Pfact;
else
    fprintf('ERROR: Bogus Method type\n');
    return;
end
Rfact              = TransPFactory();
PRfact             = GenericPRFactory(Pfact, Rfact);
PRfact.SetMaxCoarseSize(5);
Acfact             = RAPexFactory();
GSFactory          = SmootherFactory(Smoother('GaussSeidel', 2, 1));

%% This optim is not needed anymore:
%%  PRfactLvl1         = GenericPRFactory(PfactLvl1, Rfact);
%%  PRfactLvl1.SetMaxCoarseSize(100);
%%  GSFactoryLvl1      = SmootherFactory(Smoother('GaussSeidel', 2, 1, 'point'));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest, 1);
Needs.SaveAggregates = '1';
MgHierarchy.FillHierarchy(PRfact, [], Acfact, 1, numDesiredLevels);

%% This optim is not needed anymore:
%%  MgHierarchy.FillHierarchy(PRfactLvl1, [], Acfact, 1, 2);
%%  MgHierarchy.FillHierarchy(PRfact,     [], Acfact, 2, numDesiredLevels);

status = MgHierarchy.GetStatistics();

MgHierarchy.SetSmoothers(GSFactory);

%% This optim is not needed anymore:
%%  MgHierarchy.SetSmoothers(GSFactoryLvl1, 1, 1);
%%  MgHierarchy.SetSmoothers(GSFactory,     2);

NLevels = MgHierarchy.GetNumLevel();
fprintf('Levels Used = %d\n', NLevels);

oc = status.OperatorComplexity;

fprintf('OC: (2-Level) = %4.2f\n', ...
    (nnz(MgHierarchy.Levels_{1}.Get('A').GetMatrixData)  + nnz(MgHierarchy.Levels_{2}.Get('A').GetMatrixData))...
    / nnz(MgHierarchy.Levels_{1}.Get('A').GetMatrixData));

fprintf('size(level2) = %d nnz(level2) = %d\n', size(MgHierarchy.Levels_{2}.Get('A').GetMatrixData, 1), nnz(MgHierarchy.Levels_{2}.Get('A').GetMatrixData));


% HAQ
%ns = full([nullspace(:, 1);nullspace(:, 2);nullspace(:, 3);nullspace(:, 4);nullspace(:, 5);nullspace(:, 6)]);
%agg = ml('aggregate', Amat.GetMatrixData, 'PDE equations', 3, 'null space: type', 'pre-computed', ...
%       'ML output', 10, ...
%       'max levels', 2, ...
%       'aggregation: aux: enable', true, ...
%       'x-coordinates', coords(:, 1), ...
%       'y-coordinates', coords(:, 2), ...
%       'z-coordinates', coords(:, 3), ...
%             'aggregation: aux: threshold', AggDropTol, ...
%       'null space: dimension', numNS, 'null space: vectors', ns);

%keyboard;
% Compare
if(1 == 0)
    agg1 = MgHierarchy.Levels_{1}.Get('Aggregates').AggId;


    % Indices of one slice of the problem (for speed)
    %IDX = find(coords(:, 1) == 200);

    MAXX = max(coords(:, 1)); MINX = min(coords(:, 1));
    MAXY = max(coords(:, 2)); MINY = min(coords(:, 2));
    if(dim == 3), MAXZ = max(coords(:, 3)); MINZ = min(coords(:, 3));   end
    %  MATS = [0, .1, .3, .5, .7, .9, 1.1]*MAXY;
    MATS = [-.1, .33333, .66666, 1.1]*(MAXY-MINY)+MINY;

    %keyboard;

    for J = 1:NLevels-1
        if(dim == 3), coords2 = [MgHierarchy.Levels_{J}.Get('xcoords'), MgHierarchy.Levels_{J}.Get('ycoords'), MgHierarchy.Levels_{J}.Get('zcoords')];
        else coords2 = [MgHierarchy.Levels_{J}.Get('xcoords'), MgHierarchy.Levels_{J}.Get('ycoords')]; end


        % AuxMatrix after dropping at all levels (save the coarsest)
        Dmat = SAdrop(MgHierarchy.Levels_{J}.Get('AuxMatrix').GetMatrixData, 0, AggDropTol);
        figure(J+NLevels-1);
        matrix_edge_and_node_viz(Dmat, coords2);title(sprintf('Level %d connectivity', J));
        if(dim == 3), view(-90, 0);end
        hold on;
        if(dim == 3), plot3(coords2(:, 1), coords2(:, 2), coords2(:, 3), 'ro');
        else plot(coords2(:, 1), coords2(:, 2), 'ro'); end
        for I = 1:length(MATS)-1,
            if(dim == 3), patch([0, 0, 0, 0, 0], [MATS(I), MATS(I+1), MATS(I+1), MATS(I), MATS(I)], [1, 1, MAXZ+1, MAXZ+1, 1], 'k', 'FaceAlpha', 0, 'FaceColor', 'none', 'EdgeColor', 'k');
            else patch([1, 1, MAXX+1, MAXX+1, 1], [MATS(I), MATS(I+1), MATS(I+1), MATS(I), MATS(I)], [0, 0, 0, 0, 0], 'k', 'FaceAlpha', 0, 'FaceColor', 'none', 'EdgeColor', 'k');end
        end
        hold off;

        % Aggregates at all levels (save the coarsest)
        figure(J);
        agg2 = MgHierarchy.Levels_{J}.Get('Aggregates').AggId;
        if(dim == 3), plot3(coords2(:, 1), coords2(:, 2), coords2(:, 3), 'rx');view(-90, 0);
        else plot(coords2(:, 1), coords2(:, 2), 'rx'); end
        hold on;
        plot_aggregate_boxes(coords2, agg2, .4);hold off;title(sprintf('Level %d aggs', J));xlabel('x');ylabel('y');zlabel('z');
        for I = 1:length(MATS)-1,
            if(dim == 3), patch([0, 0, 0, 0, 0], [MATS(I), MATS(I+1), MATS(I+1), MATS(I), MATS(I)], [1, 1, MAXZ+1, MAXZ+1, 1], 'g', 'FaceAlpha', 0, 'FaceColor', 'none', 'EdgeColor', 'g');
            else patch([1, 1, MAXX+1, MAXX+1, 1], [MATS(I), MATS(I+1), MATS(I+1), MATS(I), MATS(I)], [0, 0, 0, 0, 0], 'g', 'FaceAlpha', 0, 'FaceColor', 'none', 'EdgeColor', 'g');end
        end
    end

    %figure(2);plot3(coords2(:, 1), coords2(:, 2), coords2(:, 3), 'rx');hold on;
    %plot_aggregate_boxes(coords2, agg2, .4);hold off;title('Level 2 aggs');xlabel('x');ylabel('y');zlabel('z');
    keyboard;
end


if(1 == 0)
    % HAQ
    ns = full([nullspace(:, 1);nullspace(:, 2);nullspace(:, 3);nullspace(:, 4);nullspace(:, 5);nullspace(:, 6)]);
    agg = ml('setup', Amat.GetMatrixData(), 'PDE equations', 3, 'null space: type', 'pre-computed', ...
        'ML output', 10, ...
        'max levels', 2, ...
        'coarse: type', 'Amesos-KLU', ...
        'print hierarchy', -1, ...
        'aggregation: aux: enable', true, ...
        'aggregation: aux: threshold', AggDropTol, ...
        'x-coordinates', coords(:, 1), ...
        'y-coordinates', coords(:, 2), ...
        'z-coordinates', coords(:, 3), ...
        'null space: dimension', numNS, 'null space: vectors', ns);

    P = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
    ml('cleanup');

    Pml = spconvert(load('Pmat_1.m'));

    Aml = spconvert(load('Amat_0.m'));

    figure(1);spy(P);title('MueMat');
    figure(2);spy(Pml);title('ML');

    quit
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful utility function from Chris' stash.
function s = cdate
v = datevec(date);
s = sprintf('%2.2d%2.2d%4d', v(2), v(3), v(1));
end