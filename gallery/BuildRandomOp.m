% Build a sparse Operator with random integer entries in the interval [a,b].
%
% syntax: [A] = BuildRandOp(n,a,b,density,seed)
%
%     n       -- number of rows
%     [a,b]   -- interval for random integers
%     density -- real number in [0,1] to control sparsity
%     seed    -- random number generator seed
%
% FIXME Right now, the diagonal values can be outside [a,b].
%
function [A] = BuildRandomOp(n,a,b,density,seed)

  if ~varexist('density'), density = 0.5; end
  if ~varexist('a'),       a = 1;         end
  if ~varexist('b'),       b = 10;        end
  if varexist('seed')
    stream = RandStream.getDefaultStream();
    reset(stream,seed);
  end
  leng = abs(b-a);
  % create matrix with random real entries.
  Adat = sprand(n,n,density);
  % put nonzeros on the diagonal
  %d = 2*abs(rand(n,1));
  %Adat = Adat + diag(d);
  lentry = Adat(n,n); Adat(n,n) = 1;
  % change all entries in matrix to integers in the desired range
  [I,J,V] = find(Adat);
  V = round(V*leng) + a;
  Adat = sparse(I,J,V);
  % Make sure there is a nonzero entry in the last row and column.
  % Otherwise the resparsifying could make the matrix nonsquare.
  if (lentry), Adat(n,n) = round(lentry*leng)+a;
  else         Adat(n,n) = 0; end
  d = diag(Adat);
  if a ~= 0, d = ~d * a;
  else       d = ~d * b;  end
  Adat = Adat + diag(d);

  newMap = Map(n, 1);

  A = Operator(Adat,newMap,newMap,@MatlabApply);
end
