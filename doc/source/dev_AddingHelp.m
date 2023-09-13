%% How to Add HTML Documentation to MueMat
%
%% Overview
% The documentation in |MueMat/doc/html| is generated with the publishing capabilities in Matlab.
%
%% Documentation source files
%
% Any .m file can be converted to an HTML file that the Matlab help browser can understand.
% For detailed information on formatting and the like, in the Matlab User Guide please see the topics
% <http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_env/f6-22451.html Publishing MATLAB Code Files> and 
% <http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_env/bruby4n-1.html Providing Your Files in the Help Browser>.
% For a summary of formatting, see <http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_env/f6-30186.html#breug1i Summary of Markup for Formatting MATLAB Files for Publishing>.
%
% Edit |example/Simple.m| to see how you go about marking up a source file. The menu 'Cell' -> 'Insert Text Markup' of the MATLAB Editor can also help you in your first steps.
%
% In the MueMat project, the main documentation is generated from the examples themselves.
% (Of course, there is also documentation associated with the source code.)
% Using examples as documentation will hopefully make it easier to keep everything up-to-date.
%
% If, however, you want to create some documentation than isn't based on code (like this file!), then the
% source should be put in the directory |MueMat/doc/source|.
%
%% Generating HTML from a single source file
%
% It's very easy to (re)generate an html page. Just type
%
%        >> makedoc 'yourFileName'
%
%% Generating HTML from all source files
%
% To generate the html pages for all sources in |MueMat/doc/source|, type
%
%        >> makealldoc
%
%% Generating just the class reference pages
%
% To generate just the class reference pages, which simply provide help for all classes in subdirectories of |MueMat/src|, type
%
%        >> makeclassdoc
%
%% Getting the HTML to show up in the Contents tab of the Matlab Help Browser
%
% To get the page to show up under MueMat toolbox in the <../images/helpwindow.png Contents> tab of the Matlab help browser,
% add the HTML filename to the file |MueMat/doc/html/helptoc.xml|.  You must restart Matlab.
%
% If you want your page to be a link on the top-level _"MueMat Toolbox"_ page, you must manually add it to
% |doc/html/muemat_product_page.html|.
%
% If you want your page to be a link on the _"Developer Resources"_ page, you must manually add it to
% |doc/source/developers.m|.
%
%% What if I don't see my changes?
%
% Try the following:
%
% * Refresh the help browser by pressing |F5|.
% * Close and reopen the help browser.
% * Restart Matlab.
