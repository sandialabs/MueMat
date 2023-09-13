%
% dblist.m  -- show some debugging context
%
% Lists a few lines of the file surrounding the line where an error is reported.
% To get more or less lines of context, change the value of n__DBG__.
% Obviously, you must be in keyboard mode.
%
% TODO It would be nice if the number of context lines you wanted could be specified
% TODO on the command line, e.g.,
% TODO     K>> dblist 5
%
dbup
n__DBG__ = 10;
[ST__DBG__,I__DBG__] = dbstack;
lmin__DBG__ = ST__DBG__(I__DBG__).line - n__DBG__;
if lmin__DBG__ < 1, lmin__DBG__ = 1; end  % starting line must be positive integer.
lmax__DBG__ = ST__DBG__(I__DBG__).line + n__DBG__;
eval([ 'dbtype ' ST__DBG__(I__DBG__).file ' ' num2str(lmin__DBG__) ':' num2str(lmax__DBG__) ]);
clear *__DBG__
