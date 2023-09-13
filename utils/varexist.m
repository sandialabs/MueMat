% Check existence of variable.
%
% Exemple: varexist('varname')
%
% The construct varexist('varname') is strictly equivalent to exist('varname','var'):
% If the variable 'varname' exist, it returns 1. 0 elsewhere.
%
% The goal of this one-liner function is to avoid the use of
% 'exist' with a single argument to test the existence of a
% variable:
% - if the variable 'varname' does not exist but a function,
% folder, or class is named 'varname', then exist('varname') return
% a value > 0 and constructs ~exist('varname') and 
% exist('varname') ~= 0 are not equivalent and error prone.
% - specifying a second argument to indicate the category to which
% the first argument belongs speed up the search.
%
% See also: exist
function [bool] = varexist(varname)
%   bool = exist(varname,'var');
    bool = evalin('caller', ['exist(''' varname ''',''var'')']);
end
