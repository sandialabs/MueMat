% Save sparse matrix into coordonnate format (I,J,V)
% Load the matrix with: M=load('M'); M=spconvert(M);
%
% Only work if there is at least one non-zero on the last row and
% columns.

function spsave(filename, M)
    [I,J,V] = find(M);
    IJV = [I,J,V];
    save(filename,'IJV','-ascii','-double');
end