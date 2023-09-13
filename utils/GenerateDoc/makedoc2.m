function makedoc2(filename, opts, outputFilename)
%
% Usage:  makedoc2 <filename> or makedoc2('<filename>',opts)
%
% Input:
%   filename (required) -- .m source file
%   opts     (optional) -- options to pass to the publish command
%
% Purpose:
%
% This utility creates an html help page from a *.m file.
% If not specified already in 'opts', the resulting html page is put
% in MueMat/doc/html.
%
% See also makealldoc, publish.
%
SetHomeDir
if ~varexist('opts') || ~isfield(opts,'outputDir')
  opts.outputDir = [MUEMAT_ROOT_DIR '/doc/html'];
end

if ~varexist('outputFilename')
  [s,outputFilename] = unix(['basename ' filename ' .m']);
  outputFilename = [outputFilename(1:end-1) 'html']; % chop()
end

txt = help2html(filename,opts);
fd = fopen([opts.outputDir '/' outputFilename],'w');
if (fd == -1), error('makedoc2', [opts.outputDir '/' outputFilename]); end
fprintf(fd,'%s',txt);
fclose(fd);
