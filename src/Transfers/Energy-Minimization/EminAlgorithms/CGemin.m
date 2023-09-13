function [P] = CGemin(A, SparsityPattern, B, P0, Nits)
%function [P] = CGemin(A, SparsityPattern, B, P0, Nits, options)
%
% This is CG applied to the problem
%
%  min  (1/2) sum (p_i)^T A p_i
%   P          i
%
%  subject to  P(i,j) = 0  if  (i,j) is not in SparsityPattern
%      and     P cnull =  fnull.
%

% Initialize
%
%blecko = (A ~= 0);
%A(options.difrows,:) = options.S(options.difrows,:) .* blecko(options.difrows,:);
%A(:,options.difrows) = options.S(:,options.difrows) .* blecko(:,options.difrows);
%A = options.S;
n       = size(A,1);
ddd     = spdiags(1 ./ diag(A),0,n,n);
%ddd(options.difrows,options.difrows) = sparse(length(options.difrows),length(options.difrows));
%scaleme = ones(n,1); scaleme(options.difrows) =  -0.00004; ddd = spdiags(scaleme,0,n,n)*ddd;
%ddd(options.difrows,options.difrows) = sparse(length(options.difrows),length(options.difrows));

% Initialization.
% Note: BBt and PInds already computed outside of the function for
%       P0 but it's not very smart to add such input parameters.
[nrows,ncols] = size(P0);
PInds =  find(reshape(SparsityPattern',[],1));
BBt = B * B';
dBBt = diag(diag(BBt));
BBt = dBBt \ BBt;
scB = dBBt \ B;
P=P0;

% CG iterations

% enforce sparsity pattern constraints
rk = -(A*P) .* SparsityPattern;
% enforce nullspace interpolation constraints
ss = Flatten(rk,PInds);
rk = BlowUp(ss-B'*(BBt \ (scB*ss) ), nrows, ncols, PInds);

%solve iteratively for the correction
% Note: we need to have P, rk, pk, and ap

fprintf('initial emin residual = %4.3e\n', full(max(max(abs(rk)))));
i = 1;
while ( (i<=Nits) && (max(max(abs(rk)))>1e-12))
  zk = ddd * rk;    % precondition by inverse diagonal
  newsum = full(sum(diag(rk'*zk)));
  if (i == 1)
     pk = zk;
  else
     betak = newsum/oldsum;
     pk = zk + pk * betak;
  end
  oldsum = newsum;

  % multiply A pk

  % enforce contraints
  ap =  (A*pk) .* SparsityPattern;
  ss=Flatten(ap,PInds);
  ap = BlowUp(ss-B'*(BBt\ (scB*ss) ),nrows,ncols,PInds);

  alphak = newsum/sum(diag(pk'*ap));

  P   = P  + alphak*pk;
  rk  = rk - alphak*ap;
  %fprintf('%d: r = %20.13e\n', i, full(max(max(abs(rk)))));
  i = i+1;
end
fprintf('final emin residual   = %4.3e\n', full(max(max(abs(rk)))));
