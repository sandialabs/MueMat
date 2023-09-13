#Usage : 
# perl AddFuncDescr.pl   ../src/Operator/Operator.m > /tmp/AddFuncDescr.m
# mv /tmp/AddFuncDescr.m ../src/Operator/Operator.m

use List::Util qw[min max];

my %hash;

$hash{'src'} = 'Object to copy'; #Object (SmootherPrototype of same type)
$hash{'mc'} = 'MATLAB Metaclass';

if ($ARGV[0] eq "MultiVector.m") {

    $hash{'numVectors'} = 'number of vectors (integer)';
    $hash{'vecLength'} = 'length of vectors (integer)';
    $hash{'ii'} = 'iith vector of the MultiVector (integer, optional, default=1:numVectors)';

    $hash{'normType'} = 'type of norm (integer, optional, default=2)';
    $hash{'vecnorm'} = 'norm(s) of vector(s) (double or array of doubles)';

    $hash{'Scalar'} = 'value (double)';
    $hash{'data'}    = 'data of the MultiVector (MATLAB matrix)';
    $hash{'udata'}   = 'data of the MultiVector (MATLAB matrix)';
    $hash{'isSparse'} = 'sparse or dense (boolean)';
    $hash{'u'} = 'a MultiVector';
    $hash{'v'} = 'a MultiVector';
    $hash{'result'} = 'a MultiVector';

    $hash{'obj'} = '';

} elsif ($ARGV[0] eq "Operator.m") {

    $hash{'label'} = 'name of Operator (string, optional, default=\'unamed\')';
    $hash{'A'} = 'a MATLAB matrix';
    $hash{'Atrans'} = 'A\' (Operator)';
    $hash{'scalar'} = 'value (double)';
    $hash{'Amat'} = 'an Operator';
    $hash{'Bmat'} = 'an Operator';
    $hash{'Cmat'} = 'an Operator';
    $hash{'Op'} = 'an Operator';
    $hash{'lambda'} = 'eigenvalue of D^{-1}.A (double)';
    $hash{'dim'} = 'get size of the dimension specified by the scalar dim (integer, optional, default=1:end)';
    $hash{'Vec'} = 'a vector (MATLAB vector)';
    $hash{'RowBlkSize'} = 'size of row block (integer)';
    $hash{'ColBlkSize'} = 'size of row block (integer)';
    $hash{'rowmap'} = 'a row map (Map)';
    $hash{'colmap'} = 'a column map (Map)';
    $hash{'diagonal'} = 'diagonal matrix (single-view Operator)';
    $hash{'data'} = 'matrix data (MATLAB matrix)';

    $hash{'view'} = 'an operator view (optional, default=defaultView)';
    $hash{'newview'} = 'a operator view (string)';
    $hash{'oldview'} = 'old operator view (string)';
    $hash{'svOp'} = 'an operator view (OperatorView)';

    $hash{'Subset'} = 'Subset of indices on which to operate (optional, default=1:end)';
    $hash{'varargin'} = 'optional parameters';

    $hash{'nrow'}  = 'number of rows (integer)';
    $hash{'ncols'} = 'number of columns (integer)';

    $hash{'OutVec'} = 'output vector';

    # TODO: NAMES NON CONSISTANT
    $hash{'applyFcn'} = 'matrix-vector multiplication function (handle)';
    $hash{'matvec'}   = 'matrix-vector multiplication function (handle)';

    $hash{'MatrixInvData'} = 'inverse matrix data (MATLAB matrix)';
    $hash{'ApplyInv'} = 'apply inverse function (handle)';

    $hash{'TorF'} = '(boolean)';
    $hash{'ViewList'} = 'list of views (cell array)';

} elsif (($ARGV[0] eq "SmootherFactory.m") || ($ARGV[0] eq "Hybrid2x2SmootherFactory.m")) {

    $hash{'Amat'} = 'matrix of the current level (Operator)';
    $hash{'PreSmootherPrototype'} = 'prototype for pre-smoothers (SmootherPrototype)';
    $hash{'PostSmootherPrototype'} = 'prototype for post-smoothers (SmootherPrototype)';
    $hash{'PreSmoo'} = 'pre-smoother (SmootherBase)';
    $hash{'PostSmoo'} = 'post-smoother (SmootherBase)';
    $hash{'Level'} = 'level of the MG hierachy (Level)';
    $hash{'Specs'} = 'specifications (CrossFactory)';

}  elsif (($ARGV[0] eq "SmootherBase.m") || ($ARGV[0] eq "SmootherBase.m")) {

    $hash{'type'} = 'identify the type of a smoother object (string)';
    $hash{'nIts'} = 'number of iterations (integer)';
    
}

if ($ARGV[0] =~ /Smoother/) {
#Smoother/
    $hash{'Amat'} = 'matrix of the current level (Operator)';
    $hash{'Level'} = 'level of the MG hierachy (Level)';
    $hash{'Specs'} = 'specifications (CrossFactory)';

    $hash{'iu'} = 'vector to smooth (could be a MultiVector or in Matlab format)';
    $hash{'irhs'} = 'right-hand side  (could be a MultiVector or in Matlab format)';
    $hash{'SolStatus'} = 'when InitGuessStatus==ALLZEROS, iu is assumed to be all zeros and initial matvec operation could be avoid (optional,default=NOTALLZEROS)';

    $hash{'Params'} = 'optional parameters for MATLAB\'s ILU (struct)';

    $hash{'NumR'} = 'A(1:NumR) == A_rr (integer)';
    $hash{'SmootherOne'} = 'smoother of the first system block (A_rr)';
    $hash{'SmootherTwo'} = 'smoother of the second system block (A_xx)';
    $hash{'MiddleTwo'}   = 'number of iterations (integer)';
    $hash{'StartTwo'}    = 'number of iterations (integer)';
    $hash{'EndTwo'}      = 'number of iterations (integer)';
    $hash{'nMainIts'}    = 'number of iterations (integer)';
    $hash{'MiddleOne'}   = 'number of iterations (integer)';

    $hash{'BlockingScheme'} = 'Block or point smoother (string=BlkDiagonal or PointDiagonal)';
    $hash{'BDiag'}          = 'diagonal';
    $hash{'string'}         = 'smoother name (string)';
    $hash{'ForwardSweep'}   = 'forward sweep (boolean)';
    $hash{'BackwardSweep'}  = 'backward sweep (boolean)';
    $hash{'Collection'}     = 'how to group blocks together for the smoothing process (Collection)';
    $hash{'Omega'}          = 'damping factor (double)';
    $hash{'nIts'}           = 'number of iteration (integer)';
    $hash{'JacobiStyle'}    = '(boolean)';

}


open (FD,$ARGV[0]);
while ($line = <FD>) {
    push(@file, $line);
}
close (FD);
$nLine = @file;

$iLine = 0;
while ($iLine <= $nLine) {
    $line = $file[$iLine];
    print $line;
    $sline = $line;
    $sline =~ s/ //g; # remove space to simplify regexp
    $sline =~ s/\t//g; # remove space to simplify regexp

    $found=0;

    # Regexp: function XXX(XXX
    if ($sline =~ /^[^%]*function(.*)\((.*)/) {

	$found = 1;
	
	# Find a line with (function [argouts] = fname(argins
	if ($sline =~ /^[^%]*function\[(.*)\]=(.*)\(([^\)]*)/) {
	    $argouts_str = $1;
	    $fname       = $2;
	    $argins_str  = $3;
	# Find a line with (function argouts = fname(argins
	} elsif ($sline =~ /^[^%]*function(.*)=(.*)\(([^\)]*)/) {
	    $argouts_str = $1;
	    $fname       = $2;
	    $argins_str  = $3;
	# Find a line with (function fname(argins
	} elsif ($sline =~ /^[^%]*function(.*)\(([^\)]*)/) {
	    $argouts_str = '';
	    $fname       = $1;
	    $argins_str  = $2;
	}
  
##
      	# if multiline arguments (argin1, argin2, ...)
        $sline = $argins_str;
        while ($sline =~ /\.\.\./) {
	    #next line is a part of function declaration
	    $iLine++;
	    $line  = $file[$iLine];
 	    print $line;
	
	    $sline = $line;
	    $sline =~ s/ //g; # remove space to simplify regexp
	    $sline =~ s/\t//g; # remove space to simplify regexp
	    
	    $argins_str .= $sline;
	}
  
        #cleanup
        $argins_str =~ s/\.\.\.//g;
        $argins_str =~ s/\n//g;
	$argins_str =~ s/\)//g;
##
        
	# Get argin list
	@argins = split(',', $argins_str);
        # remove 'this' from the list
	@argins2 = ();
	foreach my $argin (@argins) {
	    if ($argin !~ /this/) {
		push(@argins2, $argin);
	    } else {
		$addfname = 'obj.';
	    }
	}

	# Get argout list
	@argouts = split(',', $argouts_str);
        # remove 'this' from the list
	@argouts2 = ();
	foreach my $argout (@argouts) {
	    if ($argout !~ /this/) {
		push(@argouts2, $argout);
	    } else {
		#push(@argouts2, 'obj');
	    }
	}

	# Build the help line 'SYNTAX'
	$syntax = '';
	$syntax .= '[' if (@argouts2>1);
	$syntax .= $_ .", " foreach (@argouts2); if (@argouts2>0) { chop($syntax);chop($syntax); } 
	$syntax .= ']'   if (@argouts2>1);
	$syntax .= ' = ' if (@argouts2>0);
	$syntax .= $addfname . $fname;
	$syntax .= '(';
	$syntax .= $_ .", " foreach (@argins2);  if (@argins2>0)  { chop($syntax);chop($syntax); }
	$syntax .= ');';

# # # #

# how to do an union in PERL
#	my @args = (@argouts, @argins);
#	print ":\t$_\n" foreach (@args);

	# remove from argouts variables which are already in argin
	@argouts3=();
	foreach $e (@argins2) { $count{$e}++ }
	foreach my $argout (@argouts2) {
	    if (!exists $count{$argout}) {
		push(@argouts3, $argout);
	    }
	}
	@argouts2=@argouts3;
	# print ":\t$_\n" foreach (@argouts3);
        #

	@args = (@argins2, @argouts2);
# 	foreach  my $arg (@args) {
# 	    print STDERR "$arg\n";
# 	}

        # Build hash table of arg parameter list
	%hargs=();
	foreach  my $arg (@args) {
	    $hargs{$arg} = $hash{$arg};
	    if ($hash{$arg} eq '') { $unknowhash{$arg}++; }
	}

        # Description line
	$descr = ''; 	    
	$descr2 = ''; 	    
	$iLine++; # first line after function prototype

	#if next line start with %[A-Z], description already formated by this script
	# we get the description + parameter descriptions (if present) from the M-file
	if ($file[$iLine] =~ /%[A-Z]/) {
	    #print STDERR "DESCR FOUND !\n";

	    #descr
	    $descr = $file[$iLine];
	    $descr =~ s/^[ \t]*%[A-Z0-9_\(\),]*[ ]?//;
#	    if ($descr =~ /^\n$/) { $descr = ''; }

	    $iLine++; # second line after function prototype
 	    while (($file[$iLine] !~ /%   SYNTAX/) && ($file[$iLine] =~ /^[ \t]*%[ \t](.*)/)) {
		$a=$1;
		$descr .= $a."\n";
		$iLine++;
 	    }

	    #now, parse:
	    # %   SYNTAX   obj = obj.SetDinvALambda(lambda, label);
	    # %
	    # %     lambda - eigenvalue of D^{-1}.A
	    # %     label  - label

	    $iLine++; # Skip SYNTAX
	    $iLine++; # Skip blank line
	    while ($file[$iLine++] =~ /^[ \t]*%[ \t]*([^- ]*)[ \t]*-[ \t]*(.*)[ \t]*/) { #parse "varname - description"
		$a=$1;
		$b=$2;
		#print STDERR "CHECK $1\n";
		if ((exists $hargs{$a}) && ($b =~ /[A-Za-z0-9]/)) { # update description (== keep description of the m-file, not the one of the htable) only if the variable still exist on the new version of the function
		                                               # if no description in the m-file (A - ), use the description of the htable
		    $hargs{$a} = $b;
		}
	    }

	    $iLine--; #

	    #descr2
	    while ($file[$iLine++] =~ /^[ \t]*%[ \t](.*)/) {
		$descr2 .= $1."\n";
	    }
	    $iLine--; #

	} else {
	    
	    #descr
	    while ($file[$iLine++] =~ /^[ \t]*%[ \t](.*)/) {
		$descr .= $1."\n";
	    }
	    $iLine--; #
	
	}
	
#CHECKCHECKCHECK
#	print $file[$iLine];

# 	Parse description line
# 	print "--$descr"."--\n";
	if ($descr =~ /^\n$/) { $descr = ''; }

	# Build the parameter list
        # Compute the max length of a parameter name
	$lenmax = 0;
	$lenmax = max($lenmax,length($_)) foreach (keys %hargs);
	$lenmax++;

	# Build arg parameter list
	$params = '';
	foreach  my $arg (@args) {
	    $params .= $arg;
	    for ($i = length($arg); $i < $lenmax; $i++) {
		$params .= ' ';
	    }
	    $params .= '- ';
	    $params .= $hargs{$arg}."\n";
	}


	# PRINT
	$c = '      %';
	$s  = '   ';
	$s2 = '  ';

	print $c . uc($fname) . ' ';
	foreach (@list = split("\n", $descr)) {
	    print $_."\n".$c . ' ' ;
	}
	print "\n";
	if ($descr eq '') { print $c .' '."\n"; }
	#if ($descr eq '\n') { print $c .' '."\n"; }
#	else { print "BLABLA".$descr."BLAL\n"; }
#	print $descr;

	print $c . $s . 'SYNTAX   ' . $syntax ."\n";
	print $c.' ' . "\n";
	foreach (@list = split("\n", $params)) {
	    print $c . $s . $s2 . $_."\n";
	}

	if ($descr2 ne '') {
	    foreach (@list = split("\n", $descr2)) {
		print $c . " " . $_."\n";
	    }
	}

    } else {
	$iLine++;
    }
}

for my $k (keys(%unknowhash)) {
    print STDERR $k."\n";
}

