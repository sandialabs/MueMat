%% MueMat Coding Guidelines
% Coding recommendations for MueMat developers. This document presents the naming convention which are enforced in MueMat.
%
% See also : <http://www.mathworks.com/matlabcentral/fileexchange/2529 MATLAB Programming Style Guidelines by Richard Johnson>
% 
%% Classes
%  - Class names should be nouns in CamelCase with the first letter capitalized.
%
%      Example(s): SmootherFactory
%
%  - Class methods should be verbs in CamelCase with the first letter capitalized.
%
%      Example(s): GetProperty(), SetProperty(), Build().
%
%  - Class properties ended with an underscore (_)
%
%  - Constants should be all uppercase (currently, we don't use underscore to separate words in mue_include but we could).
%
%% Variables and parameters
% 
% - Variables should be named in camelCase with the first letter lowercase.
%
%      Exceptions:
%
%        * Name of matrices are upper case.
%
%           Example(s): P / A / Pdata / APpattern / APnode / Mdata / ...
%
%        * Name of variables storing instances of factories are also upper case.
%
%           Example(s): SmootherFact
%        
% - The prefix 'n' should be used for variables representing the number of objects.
%
%      Example(s): nRows
%
% - Methods related to variables representing a number should use "Num" on their names.
%
%      Example(s): SetNumIterations(nIts)
%
% - Convention on pluralization: 's' suffix.
%
%      Example(s): nRows, nIts
%
% - Convention on acronyms: first letter capitalized, all the rest lowercase.
%
%      Example(s): Nnz, Dof, Rhs
% 
% - Convention on compound: no capitalization of the first letter of second stem.
%
%      Example(s): nullspace, multigrid    
%
% - Please choose meaningful names. We also like names that reflect variable function. 
%
%      Example(s): 
%
%       * CoarseDofPerNode (vs. nCDofs)
%      
%       * finestLevel (vs. finest)
%      
%       * aggregates (vs. aggInfo or aggWidget)
%
% - List of allowed abbreviations:
%
%        * blk for block
%
%        * const for constant
%
%        * sol for solution (sol = mg.Iterate())
%
%        * its for iterations
%
%        * null for nullspace
%
%% Indentation, spacing, tabs
%
% Please do *not* insert tab characters into the code. Instead, use the equivalent number of spaces.
% Here are instructions for causing the tab key to insert spaces rather than a tab character:
%
% *vim*   --   add the line |set expandtab| to your |.vimrc| file.
%
% *emacs* --   add the line |(setq-default indent-tabs-mode nil)| to your |.emacs| file.
%
% The emacs hint comes from this <http://www.gnu.org/software/emacs/emacs-lisp-intro/html_node/Indent-Tabs-Mode.html document>.
