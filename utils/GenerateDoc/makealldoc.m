function makealldoc(opts)
%
% Usage:  makealldoc or makealldoc(opts)
%
% Input:
%   opts     (optional) -- options to pass to the publish command
%
% Purpose:
%
% This utility generates all the html help page from .m files.
% TODO: for the moment, this script generates help from *.m files of doc/source *only*
%
% If not specified already in 'opts', the resulting html page is put
% in MueMat/doc/html.
%
% See also makedoc, publish.
%
SetHomeDir
if ~varexist('opts') || ~isfield(opts,'outputDir')
  opts.outputDir = [MUEMAT_ROOT_DIR '/doc/html'];
end

oldpath=path;
path([path ':' MUEMAT_ROOT_DIR '/doc/source/']);

filelist = dir([MUEMAT_ROOT_DIR '/doc/source/' '*.m']);
for i=1:size(filelist,1)
    publish(['/doc/source/' filelist(i).name],opts);
end

path(oldpath);