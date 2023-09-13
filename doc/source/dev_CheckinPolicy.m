%% Policy/Process for Safe Checkin Testing
%
%% Non-Regression testing
% In order to maintain the stability of MueMat, the shell script 'non-regression.sh' should always be used to test modifications before any checkin that changes source code. 
% This script executes a batch of tests. Outputs of tests are compared against outputs of previous MueMat version.
%
%% Setting environment
%
% The shell script can be found in 'test/NonRegression/'. 
%
% PreRequirement:
%
% - GNU/parallel is used to run test cases in parallel (by using one instance of MATLAB per processor).
%   GNU/parallel can be downloaded at: http://www.gnu.org/software/parallel/
%
%  Quick install guide:
%    tar -xjf parallel-*.tar.bz2
%    cd parallel-*
%    ./configure
%    make -j; sudo make install
%
% - tkdiff is used by default to present results of the non-regression testing.
%
%% Script usage
%
%   ./non-regression.sh: run tests and compare results 
%                        against the last version of output 
%                        available in 'REF/'.
%
% A tkdiff window is open for each test output that differ from the previous version. You can check visually if differences are significative.
%
%% Managing version of reference
%
% The directory 'current-version/' is temporarily used to store outputs of your test runs. This directory can be removed after non-regression testing.
% The directory 'REF/' stores outputs of previous MueMat versions.
%
% If you fix a bug and believe that changes on the results of tests are OK, you can commit a new version of output reference in REF/ by using the following command:
%
%   mv current-version REF/`date +%F`-`git rev-parse HEAD | tail -c 10`
%
%% Adding new tests
%
% The non-regression script executes the MATLAB commands defined in 'config.txt'.
%
%% Advanced usages
%
% - non-regression.sh uses internally two sub-scripts that can be run independently: 
%
%   1) non-regression-run.sh:   run tests.
%   2) non-regression-check.sh: compare results against reference output.
%
% - Advanced options are listed on the header of shell scripts.