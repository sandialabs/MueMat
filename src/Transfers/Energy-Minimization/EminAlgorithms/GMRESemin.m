function [P] = GMRESemin(A, SparsityPattern, B, P0, Nits)
%gmres_r   GMRES iteration with right preconditioning
%   [x] = gmres_r(mparams, b, params, x0)
%   input
%          amat
%          mparams      structure defining preconditioning matrix
%          b            right-hand side vector
%          params       three-dimensional vector to control iteration
%            params(1) = relative residual reduction factor
%            params(2) = max number of iterations
%            params(3) (Optional) = reorthogonalization method
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%          x0           initial iterate
%   output
%          x            computed (approximate) solution
%   IFISS function: HCE; 15 March 2005.


%Modification of gmres.m written by Tim Kelley.
%Copyright 1994 C. T. Kelley.
%Reproduced and distributed with permission of the copyright holder.
%Modifications to handle IFISS data structures identified in line.

% GMRES linear equation solver
% Implementation following Saad-Schultz
%
% C. T. Kelley, July 10, 1994
%
% Requires givapp.m
% global bcs;


%
% initialization
%
n = size(A,1);
%[nfine,ncoarse] = size(SparsityPattern);
%ntotal = nfine*ncoarse;
%nulldim = size(NSfine,2);
%ddd     = spdiags(1 ./ (abs(A)*ones(n,1)),[0],n,n);
ddd     = spdiags(1 ./ (ones(n,1)),0,n,n);
%ddd     = spdiags(1 ./ diag(A),[0],n,n);
% if (nulldim ~= 0)
%    NSCoarse    = NSfine(roots,:);
% else
%    NSCoarse    = [];
% end

reorth=0;

% Set the intial P
% Initialization.
% Note: BBt and PInds already computed outside of the function for
%       P0 but it's not very smart to add such input parameters.
[nrows,ncols] = size(P0);
PInds =  find(reshape(SparsityPattern',[],1));
BBt = B * B';
dBBt = diag(diag(BBt));
BBt = dBBt \ BBt;
P=P0;


%
%
h=zeros(Nits);
c=zeros(Nits+1,1);
s=zeros(Nits+1,1);
r = -ddd*A*P;
r =  r .* SparsityPattern;

%fprintf('taking out constraint stuff and just zeroing out cpt update!!!\n');
%fprintf('putting in crazy bc stuff\n');
%if (~isempty(bcs) ) old = r(bcs,:); end
ss = Flatten(r,PInds); r = BlowUp(ss-B'*(BBt\(dBBt \(B*ss))),nrows,ncols,PInds);
%if (~isempty(bcs) ) r(bcs,:) = old; end
%r = NoCptUpdate(r, SparsityPattern, NSCoarse, BtBFact, nPDE,roots);
%fprintf('QA Gmres:\n');
rho=norm(r,'fro');         %fprintf('  %5i    %20.13e \n',0,rho);
g=rho*eye(Nits+1,1);
%
% test for termination on entry
%
%
v.num(1).matrix = r/rho;
%beta=rho;
k=0;
%
% GMRES iteration
%
while(k < Nits)
    k=k+1;
    v.num(k+1).matrix = ddd*A*v.num(k).matrix;
    v.num(k+1).matrix =  v.num(k+1).matrix .* SparsityPattern;
%if (~isempty(bcs) ) old = v.num(k+1).matrix(bcs,:); end
ss = Flatten(v.num(k+1).matrix,PInds); v.num(k+1).matrix=BlowUp(ss-B'*(BBt\ (dBBt \(B*ss))),nrows,ncols,PInds);
%if (~isempty(bcs) ) v.num(k+1).matrix(bcs,:) = old; end
%    v.num(k+1).matrix = NoCptUpdate(v.num(k+1).matrix,SparsityPattern,...
%                               NSCoarse, BtBFact, nPDE,roots);
    normav=norm(v.num(k+1).matrix,'fro');
%
% Modified Gram-Schmidt
%
    for j=1:k
       h(j,k)=sum(diag(v.num(j).matrix'*v.num(k+1).matrix));% inner product of j
                                                            % and k+1 vecs

       v.num(k+1).matrix = v.num(k+1).matrix - h(j,k)*v.num(j).matrix;
    end
    h(k+1,k)= norm(v.num(k+1).matrix,'fro');
    normav2=h(k+1,k);
%
% Reorthogonalize?
%
if  (reorth == -1 && normav + .001*normav2 == normav)
    for j=1:k
        hr=v(:,j)'*v(:,k+1);                 % hr=v(:,k+1)'*v(:,j);
        h(j,k)=h(j,k)+hr;
        v(:,k+1)=v(:,k+1)-hr*v(:,j);
    end
    h(k+1,k)=norm(v(:,k+1));
end
%
%   watch out for happy breakdown
%
    if (h(k+1,k) ~= 0)
        v.num(k+1).matrix = v.num(k+1).matrix/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation
%
    if k > 1
        h(1:k,k)=givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
    nu=norm(h(k:k+1,k));
    if nu~=0
        c(k)=conj(h(k,k)/nu);             %   c(k)=h(k,k)/nu;  % Change 6/3/97
        s(k)=-h(k+1,k)/nu;
        h(k,k)=c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k)=0;
        g(k:k+1)=givapp(c(k),s(k),g(k:k+1),1);
    end
%
% Update residual norm
%
    %rho=abs(g(k+1));         %fprintf('  %5i     %8.4e \n',k,rho);
end
%

y=h(1:k,1:k)\g(1:k);
for i=1:k
   P = P + v.num(i).matrix*y(i);
end

end

function vrot=givapp(c,s,vin,k)
  %givapp   apply a sequence of Givens rotations
  %   input
  %          c, s      vectors of length k-1 defining rotations
  %          vin       vector of length k to which rotations are applied
  %          k         k-1 = number of rotations
  %   output
  %          tranformed vector after rotations are applied
  % called by gmres_r

  %
  %  C. T. Kelley, July 10, 1994
  %Copyright 1994 C. T. Kelley.
  %Reproduced and distributed with permission of the copyright holder.
  %
  % This code comes with no guarantee or warranty of any kind.
  %
  %  function vrot=givapp(c, s, vin, k)
  %
  vrot=vin;
  for i=1:k
      w1=c(i)*vrot(i)-s(i)*vrot(i+1);        % Change on next line, 6/3/97
      w2=s(i)*vrot(i)+conj(c(i))*vrot(i+1);  % w2=s(i)*vrot(i)+c(i)*vrot(i+1);
      vrot(i:i+1)=[w1,w2];
  end
end %givapp()
