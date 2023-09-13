% function [K, M, Nullspace] = BuildHelmholtz2D(waves,npts,medium)
%
% Builds the stiffness matrix, mass matrix, and near-nullspace for the 2D Helmholtz equation
% with PML on the unit square [0,1] x [0,1]
%
% Parameters:
% 'waves'  - Number of wavelengths across the interval.
% 'npts'   - Points per wavelength.
% 'medium' - type of material
%
% Output:
% K         - stiffness matrix
% M         - mass matrix
% Nullspace - near-nullspace of Helmholtz operator at given frequency

function [K, M, Nullspace, Nx, Ny] = BuildHelmholtz2D(waves,npts,medium)

  N=waves*npts;
  h=1/(N+1);
  gs=[0:h:1]';
  [xx,yy] = ndgrid(gs,gs);
  if(medium==0)
    % homogeneous medium
    c = ones(size(xx));
  elseif(medium==1)
    % converging lens
    c = ones(size(xx)) - 0.5*exp(-32*((xx-1/2).^2+(yy-1/2).^2));
  elseif(medium==2)
    % waveguide
    c = ones(size(xx)) - 0.5*exp(-32*((xx-1/2).^2));
  elseif(medium==3)
    % random medium
    c = 1+0.4*rand(8,8);
    c = interpft(interpft(c,numel(gs),1),numel(gs),2);
  else
    % homogeneous medium
    c = ones(size(xx));
  end
  % normalize coefficients around 1
  cmid = (min(c(:))+max(c(:)))/2;  c = c/cmid;

  lambda = 1/waves;
  LB = lambda;
  RB = 1-lambda;
  delta = 2.0;
  sig1 = zeros(size(gs));
  gd = find(gs<LB);  sig1(gd) = delta*abs((gs(gd)-LB)/lambda).^2;
  gd = find(gs>RB);  sig1(gd) = delta*abs((gs(gd)-RB)/lambda).^2;
  sig2 = zeros(size(gs));
  gd = find(gs<LB);  sig2(gd) = delta*abs((gs(gd)-LB)/lambda).^2;
  gd = find(gs>RB);  sig2(gd) = delta*abs((gs(gd)-RB)/lambda).^2;
  s1 = 1+i*sig1;
  s2 = 1+i*sig2;

  % setup matrices
  [K,M]=setupHelmholtz2D(h,c,s1,s2);

  % compute near-nullspace - plane waves
  omega=2*pi/lambda;
  xx=xx(2:end-1,2:end-1);
  yy=yy(2:end-1,2:end-1);
  Nullspace(:,1)=exp(i*omega*xx(:));
  Nullspace(:,2)=exp(-i*omega*xx(:));
  Nullspace(:,3)=exp(i*omega*yy(:));
  Nullspace(:,4)=exp(-i*omega*yy(:));
  Nx=size(xx,1);
  Ny=size(yy,1);

end

function [K,M] = setupHelmholtz2D(h,c,s1,s2)

  h2 = (h*h);
  [P1,P2] = size(c);
  idx = zeros(P1,P2);
  idx(2:P1-1,2:P2-1) = reshape(1:(P1-2)*(P2-2), P1-2, P2-2);
  
  MD1 = 2:P1-1;
  LF1 = 1:P1-2;
  RT1 = 3:P1;
  
  MD2 = 2:P2-1;
  LF2 = 1:P2-2;
  RT2 = 3:P2;
  
  s1 = s1(:)*ones(1,P2);
  s2 = ones(P1,1)*transpose(s2(:));
  
  Il = idx(MD1,MD2);
  Jl = idx(LF1,MD2);
  Sl = -(1/h2)*(s2(MD1,MD2)./s1(MD1,MD2)+s2(LF1,MD2)./s1(LF1,MD2))/2;
  
  Ir = idx(MD1,MD2);
  Jr = idx(RT1,MD2);
  Sr = -(1/h2)*(s2(MD1,MD2)./s1(MD1,MD2)+s2(RT1,MD2)./s1(RT1,MD2))/2;
  
  Iu = idx(MD1,MD2);
  Ju = idx(MD1,LF2);
  Su = -(1/h2)*(s1(MD1,MD2)./s2(MD1,MD2)+s1(MD1,LF2)./s2(MD1,LF2))/2;
  
  Id = idx(MD1,MD2);
  Jd = idx(MD1,RT2);
  Sd = -(1/h2)*(s1(MD1,MD2)./s2(MD1,MD2)+s1(MD1,RT2)./s2(MD1,RT2))/2;
  
  Ik = idx(MD1,MD2);
  Jk = idx(MD1,MD2);
  Sk = -(Sl+Sr+Su+Sd);
  Sm = (s1(MD1,MD2).*s2(MD1,MD2)).*(c(MD1,MD2).^2);

  Is = [Il(:); Ir(:); Id(:); Iu(:); Ik(:)];
  Js = [Jl(:); Jr(:); Jd(:); Ju(:); Jk(:)];
  Ss = [Sl(:); Sr(:); Sd(:); Su(:); Sk(:)];
  
  gd = find(Js>0);
  Is = Is(gd);
  Js = Js(gd);
  Ss = Ss(gd);
  K  = sparse(Is(:),Js(:),Ss(:));
  M  = sparse(Ik(:),Jk(:),Sm(:));

end
