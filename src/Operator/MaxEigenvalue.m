function lambda = MaxEigenvalue(Afun, n, varargin)

  opts.tol=1e-1; opts.maxit = 20; opts.p=10; opts.disp = 0; opts.isreal=false;
  opts.v0=rand(RandStream.create('mrg32k3a','NumStreams',1),n,1);
  if opts.p > n, opts.p = n; end

  lambda = abs(eigs(Afun,n,1,'LM',opts, varargin{:}));

  % keep only 2 digits of accuracy to improve chances of
  % reproducability between repeated runs.
  lambda = 1.05*ceil(100.*lambda)/100;
end
