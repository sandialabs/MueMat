% function [K, M, Nullspace] = BuildHelmholtz3D(waves,npts,medium)
%
% Builds the stiffness matrix, mass matrix, and near-nullspace for the 3D Helmholtz equation
% with PML on the unit cube [0,1] x [0,1] x [0,1]
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

function [K, M, Nullspace, Nx, Ny, Nz] = BuildHelmholtz3D(waves,npts,medium)

  N=waves*npts;
  h=1/(N+1);
  gs=[0:h:1]';
  [xx,yy,zz] = ndgrid(gs,gs,gs);
  if(medium==0)
    % homogeneous medium
    c = ones(size(xx));
  elseif(medium==1)
    % converging lens
    c = ones(size(xx)) - 0.4*exp(-32*((xx-1/2).^2+(yy-1/2).^2+(zz-1/2).^2));
  elseif(medium==2)
    % waveguide
    c = ones(size(xx)) - 0.4*exp(-32*((xx-1/2).^2+(yy-1/2).^2));
  elseif(medium==3)
    % random medium
    c = 1+0.4*rand(8,8,8);
    c = myinterpft3(c,numel(gs));
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
  sig3 = zeros(size(gs));
  gd = find(gs<LB);  sig3(gd) = delta*abs((gs(gd)-LB)/lambda).^2;
  gd = find(gs>RB);  sig3(gd) = delta*abs((gs(gd)-RB)/lambda).^2;
  s1 = 1+i*sig1;
  s2 = 1+i*sig2;
  s3 = 1+i*sig3;

  % setup matrices
  [K,M]=setupHelmholtz3D(h,c,s1,s2,s3);

  % compute near-nullspace - plane waves
  omega=2*pi/lambda;
  xx=xx(2:end-1,2:end-1);
  yy=yy(2:end-1,2:end-1);
  zz=zz(2:end-1,2:end-1);
  Nullspace(:,1)=exp(i*omega*xx(:));
  Nullspace(:,2)=exp(-i*omega*xx(:));
  Nullspace(:,3)=exp(i*omega*yy(:));
  Nullspace(:,4)=exp(-i*omega*yy(:));
  Nullspace(:,5)=exp(i*omega*zz(:));
  Nullspace(:,6)=exp(-i*omega*zz(:));
  Nx=size(xx,1);
  Ny=size(yy,1);
  Nz=size(zz,1);

end

function [K,M] = setupHelmholtz3D(h,c,s1,s2,s3)

  h2=h*h;
  [P1,P2,P3] = size(c);
  idx = zeros(P1,P2,P3);
  idx(2:P1-1,2:P2-1,2:P3-1) = reshape(1:(P1-2)*(P2-2)*(P3-2), [P1-2, P2-2, P3-2]);
  
  MD1 = 2:P1-1;
  LF1 = 1:P1-2;
  RT1 = 3:P1;
  
  MD2 = 2:P2-1;
  LF2 = 1:P2-2;
  RT2 = 3:P2;
  
  MD3 = 2:P3-1;
  LF3 = 1:P3-2;
  RT3 = 3:P3;
  
  s1 = reshape(s1,[P1,1,1]);  s1 = repmat(s1,[1,P2,P3]);
  s2 = reshape(s2,[1,P2,1]);  s2 = repmat(s2,[P1,1,P3]);
  s3 = reshape(s3,[1,1,P3]);  s3 = repmat(s3,[P1,P2,1]);
  
  Il = idx(MD1,MD2,MD3);
  Jl = idx(LF1,MD2,MD3);
  Sl = -(1/h2)*( s2(MD1,MD2,MD3).*s3(MD1,MD2,MD3)./s1(MD1,MD2,MD3) + s2(LF1,MD2,MD3).*s3(LF1,MD2,MD3)./s1(LF1,MD2,MD3) )/2;
  
  Ir = idx(MD1,MD2,MD3);
  Jr = idx(RT1,MD2,MD3);
  Sr = -(1/h2)*( s2(MD1,MD2,MD3).*s3(MD1,MD2,MD3)./s1(MD1,MD2,MD3) + s2(RT1,MD2,MD3).*s3(RT1,MD2,MD3)./s1(RT1,MD2,MD3) )/2;
  
  Iu = idx(MD1,MD2,MD3);
  Ju = idx(MD1,LF2,MD3);
  Su = -(1/h2)*( s1(MD1,MD2,MD3).*s3(MD1,MD2,MD3)./s2(MD1,MD2,MD3) + s1(MD1,LF2,MD3).*s3(MD1,LF2,MD3)./s2(MD1,LF2,MD3) )/2;
  
  Id = idx(MD1,MD2,MD3);
  Jd = idx(MD1,RT2,MD3);
  Sd = -(1/h2)*( s1(MD1,MD2,MD3).*s3(MD1,MD2,MD3)./s2(MD1,MD2,MD3) + s1(MD1,RT2,MD3).*s3(MD1,RT2,MD3)./s2(MD1,RT2,MD3) )/2;
  
  If = idx(MD1,MD2,MD3);
  Jf = idx(MD1,MD2,LF3);
  Sf = -(1/h2)*( s1(MD1,MD2,MD3).*s2(MD1,MD2,MD3)./s3(MD1,MD2,MD3) + s1(MD1,MD2,LF3).*s2(MD1,MD2,LF3)./s3(MD1,MD2,LF3) )/2;
  
  Ib = idx(MD1,MD2,MD3);
  Jb = idx(MD1,MD2,RT3);
  Sb = -(1/h2)*( s1(MD1,MD2,MD3).*s2(MD1,MD2,MD3)./s3(MD1,MD2,MD3) + s1(MD1,MD2,RT3).*s2(MD1,MD2,RT3)./s3(MD1,MD2,RT3) )/2;
  
  Ik = idx(MD1,MD2,MD3);
  Jk = idx(MD1,MD2,MD3);
  Sk = -(Sl+Sr+Su+Sd+Sf+Sb);
  Sm = (s1(MD1,MD2,MD3).*s2(MD1,MD2,MD3).*s3(MD1,MD2,MD3)).*(c(MD1,MD2,MD3).^2);
  
  Is = [Il(:); Ir(:); Id(:); Iu(:); If(:); Ib(:); Ik(:)];
  Js = [Jl(:); Jr(:); Jd(:); Ju(:); Jf(:); Jb(:); Jk(:)];
  Ss = [Sl(:); Sr(:); Sd(:); Su(:); Sf(:); Sb(:); Sk(:)];
  
  gd = find(Js>0);
  Is = Is(gd);
  Js = Js(gd);
  Ss = Ss(gd);
  K  = sparse(Is(:),Js(:),Ss(:));
  M  = sparse(Ik(:),Jk(:),Sm(:));

end
