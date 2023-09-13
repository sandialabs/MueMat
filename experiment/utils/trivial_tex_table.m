% function trvial_tex_table(fname,prf,DATA,DFOR,ROW_LBL,COL_LBL,CAPTION)
% Generates a simple table in LaTeX w/o blocking
%
% Input:
% fname      - output file name
% prf        - prefix for the item label
% DATA       - 2D array of data, in (ROW,COLUMN) format
% DFORMAT    - printf format string for output data [array or scalar]
% ROW_LBL    - Labels for the rows, ROW_LBL{1} should be the 
%              column label for the rows themselves.
% COL_LBL    - Labels for the columns.
% CAPTION    - The caption for the figure [string]
%
% [AST_MAX]    - Put asterisks next to data entries that equal
%                or exceed AST_MAX. [default=inf].
% [AST_CAPTION]- Label to add to the caption if AST_MAX is equalled
%                or exceeded. [default=not used].
% by: Chris Siefert <csiefer@sandia.gov>
% Version History:
% 0.1 10/16/2006 - Initial Version <csiefer@sandia.gov>
% 0.2 02/06/2008 - Added AST_MAX / AST_CAPTION <csiefer@sandia.gov>
% 0.3 07/28/2009 - Adding item label.

function trivial_tex_table(fname,prf,DATA,DFORMAT,ROW_LBL,COL_LBL,CAPTION,varargin)
% Optional Arguments
if(length(varargin)>=2) AST_MAX=varargin{1};AST_CAPTION=varargin{2};
else AST_MAX=inf;end
AST_TRIGGERED=0;
  
% Calculate DFOR
if(length(DFORMAT)==1), DFOR=repmat(DFORMAT,1,length(COL_LBL));
else DFOR=DFORMAT;end
if(length(DFOR)~=length(COL_LBL)),fprintf('Error: DFORMAT / COL_LBL size mismatch (%d,%d)!\n',length(DFOR),length(COL_LBL));return;end

% Output File Header
f=fopen(fname,'w');
NCOLS=size(DATA,2);
PREF=sprintf('%s',prf);
if NCOLS>6
    fprintf(f,'\\begin{table}\\begin{center}\\resizebox{\\columnwidth}{!}{\\begin{tabular}{|r||');
else
    fprintf(f,'\\begin{table}\\begin{center}\\small{\\begin{tabular}{|r||');
end

% Formatting data w/ vinterlace
for I=1:length(COL_LBL),
  fprintf(f,'r|');
end  
fprintf(f,'}\n');
fprintf(f,'\\hline\n ');

% Print the column Labels
fprintf(f,'\\multicolumn{1}{|c||}{%s} ',ROW_LBL{1}); ROW_LBL=ROW_LBL(2:length(ROW_LBL));
for I=1:length(COL_LBL),fprintf(f,'& \\multicolumn{1}{|c|}{%s} ',COL_LBL{I});end
fprintf(f,'\\\\\n\\hline\\hline\n');  

% Print the data
for J=1:length(ROW_LBL),    
  fprintf(f,'%s ',ROW_LBL{J});
  for K=1:length(COL_LBL),
    if(iscell(DATA(J,K))), fprintf(f,sprintf('& %s',DFOR{K}),DATA{J,K});
    elseif(DATA(J,K)==0), fprintf(f,'& ---');      
    elseif(DATA(J,K)>=AST_MAX), fprintf(f,'& $\\ast$');AST_TRIGGERED=1;
    else fprintf(f,sprintf('& %s',DFOR{K}),DATA(J,K));end  
  end
  fprintf(f,' \\\\\\hline\n');
end
%fprintf(f,'\\hline\n');
fprintf(f,'\\end{tabular}}\\end{center}\n');

if(AST_TRIGGERED) fprintf(f,'\\caption{%s  %s}\n',CAPTION,AST_CAPTION);
else fprintf(f,'\\caption{%s}\n',CAPTION);end
fprintf(f,'\\label{tbl:%s}\n',prf);
fprintf(f,'\\end{table}\n');
fclose(f);

