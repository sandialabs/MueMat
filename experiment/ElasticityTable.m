% Generates charts for Elasticity Test

% Usage:
% - ElasticityTable()           : generates charts for each etest*.mat of the
%                                 current directory
% - ElasticityTable('05302012') : generates charts for etest*05302012*.mat
% - ElasticityTable('2d')       : generates charts only for 2d runs
% - ElasticityTable('2d*s10.')  : generates charts only for dim=2, stretch=10

% Output files: etest*.tex, results.txt

% FIXME: LaTeX files are overwritten if same dim/stretch with different dates.

function ElasticityTable(PATTERN, EXPLIST, COL_LBL)
if ~varexist('PATTERN'), PATTERN = ''; end

%bascic
%EXPLIST={'sa-norotate', 'sa', 'emin', 'emin-sa'};
%COL_LBL={'SA-NR', 'SA', 'EMIN', 'EMIN-SA'};
% uncomment to print more columns
%EXPLIST={'ml-sa-nr', 'ml-sa', 'emin', 'sa-norotate', 'sa'};
%COL_LBL={'SA-NR', 'SA', 'EMIN', 'SA-NR-MueMat', 'SA-MueMat'};

%EXPLIST={'sa-norotate', 'sa', 'emin', 'emin-sa'};
%COL_LBL={'SA-NR', 'SA', 'EMIN', 'EMIN-SA'};

% A first selection of columns
%EXPLIST={'nonsmoo-nr', 'sa-norotate', 'sa-nr-nf', 'emin-nr', 'emin', 'nonsmoo', 'sa', 'sa-nf', 'emin-sa'};
%COL_LBL={'NSA-NR', 'SA-NR', 'SA-NF-NR', 'EMIN-NR', 'EMIN-??', 'NSA', 'SA', 'SA-NF', 'EMIN'};

%
%
%

% For stretch>0
if ~varexist('EXPLIST'), EXPLIST={'sa-norotate', 'sa-nr-nf', 'emin-nr', 'nonsmoo', 'sa', 'sa-nf', 'emin-sa'}; end
if ~varexist('COL_LBL'), COL_LBL={'SA-NR', 'SA-NF-NR', 'EMIN-NR', 'NSA', 'SA', 'SA-NF', 'EMIN'}; end

% For stretch=0
%EXPLIST={'sa-norotate', 'emin-nr', 'nonsmoo', 'sa', 'emin-sa'};
%COL_LBL={'SA-NR', 'EMIN-NR', 'NSA', 'SA', 'EMIN'};

fd = fopen('results.txt', 'w');

list = dir(['etest*' PATTERN '*.mat']);
for i = 1:length(list)
    filename = list(i).name;

    s = regexp(filename, '\.', 'split');   date = s{length(s)-1};
    s = regexp(filename, 't|d', 'split');  dim = str2num(s{3});
    s = regexp(filename, 's|\.', 'split'); if strcmp(s{4}, 'mat') stretch = 1; else stretch = str2num(s{4}); end

    ElasticityTxtTable(dim, stretch, date);
    ElasticityTxtTable(dim, stretch, date, fd);
    ElasticityTexTable(dim, stretch, date, EXPLIST, COL_LBL, 100);
end

fclose(fd);

end