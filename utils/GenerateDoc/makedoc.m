function makedoc(filename,opts)
%
% Usage:  makedoc <filename> or makedoc('<filename>',opts)
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
publish(filename,opts);
