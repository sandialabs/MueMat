% LaTeX charts for Elasticity Test

% EXPLIST = list of columns to display
% COL_LBL = columns label

function ElasticityTable(dim,stretch,DATE,EXPLIST,COL_LBL,maxIt,tag)

if ~varexist('dim'),     dim=2; fprintf('DIM=default=%d\n', dim); end
if dim ~= 2 && dim ~= 3, fprintf('Invalid dimension\n'); keyboard; end
if ~varexist('stretch'), stretch=1; end % default = no stretching
if ~varexist('EXPLIST'), EXPLIST={'ml-sa-nr', 'ml-sa', 'emin'}; COL_LBL={'SA-NR', 'SA', 'EMIN'}; end
if ~varexist('COL_LBL'), COL_LBL=EXPLIST; end
if ~varexist('tag'), tag=[]; end

% Load File: (ITERS, LABEL, NLEVELS, NLIST, OC, RESID)
if(stretch==1), data = load(sprintf('etest%dd-filter3_304.%s.mat',dim,DATE));
else data = load(sprintf('etest%dd-filter3_304.s%d.%s.mat',dim,stretch,DATE));end

% Extract columns that are displayed (EXPLIST parameter)
% /!\ only done for ITERS and OC
ITERS = -ones(size(data.NLIST,2),size(EXPLIST,2));
OC = -ones(size(data.NLIST,2),size(EXPLIST,2));
for I=1:size(EXPLIST,2),
    idx = ismember(data.LABEL, EXPLIST{I})==1;
    if (sum(idx) == 1)
        ITERS(:,I) = data.ITERS(:,idx);
        OC(:,I) = data.OC(:,idx);
    end
end

% Labels
FPREF=sprintf('etest%dd_s%d%s.tex',dim,stretch,tag);
CLBL=sprintf('tbl:etest%dd_s%d%s',dim,stretch,tag);

% First column
ROW_LBL{1}='Size';
for I=1:size(ITERS,1),
    % Problem size
    ROW_LBL{I+1}=sprintf('$%3d^%d$',data.NLIST(I),dim);
end

% Iterations & OC
for I=1:size(ITERS,1),
    for J=1:size(ITERS,2),
        if (ITERS(I,J) >= maxIt)
            DATA{I,J}=sprintf('-- (%4.2f) ',OC(I,J));
        else
            DATA{I,J}=sprintf('%2d (%4.2f) ',ITERS(I,J),OC(I,J));
        end
    end
end

% Label type
if(dim==2), Dstring='2D plane stress';
else Dstring='3D';end

if(stretch==1), Tstring='isotropic elasticity';
else Tstring=sprintf(['anisotropic elasticity with %d:1 mesh stretching in ' ...
                     'the $x$ direction'],stretch);end

CAPTION=sprintf(['Iteration counts and operator complexities (in parenthesis) ' ...
                 'of smoothed agregation w/o rotational modes (SA-NR), ' ...
                 'smoothed aggregation (SA) and energy minimization (EMIN) ' ...          
                 'for %s %s.'],Dstring,Tstring);

DFORMAT = cell(1,size(ITERS,2));             
for i = 1:size(ITERS,2), DFORMAT{i} = '%s'; end;

trivial_tex_table(FPREF,CLBL,DATA,DFORMAT,ROW_LBL,COL_LBL,CAPTION);
