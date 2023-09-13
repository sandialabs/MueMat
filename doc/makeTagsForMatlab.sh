#!/bin/sh
# To make tags for Matlab files.  Chris S. found this hint at http://www.hausmilbe.net.
ctags -R --langdef=matlab --langmap=matlab:.m --regex-matlab="/^function([] A-Za-z0-9,_[]+=[ ]?|[ ]+)([^.(]+).*$/\2/f,function/"
