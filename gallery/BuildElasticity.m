% function [Amat, NullSpace, NODES] = BuildElasticity(n, E, nu,[stretch])
% Build a matrix (and near-nullspace) for 23D linear elasticity.
%
% This uses a tri-linear elements on a regular mesh (allows for
% mesh stretching, if requested)
%
% Parameters:
% 'NELS' - Number of elements in each dimension
% 'E'    - Elastic modulus.
% 'nu'   - Poisson's ratio.
% Optional Parameters:
% 'stretch' - A 2/3-vector containing the mesh streching in each
%             dimension ([1,1,1] = isotropic mesh).
%
% Example:
% [M,nullspace]=BuildElasticity(2 ,1e5, 0.3);
%
% References:
% [1] R.D. Cook, D.S. Malkus and M.E. Plesha
%     "Concepts and Applications of Finite Element Analysis, 3rd ed.
%
% See also: BuildElasticity2D

function [Amat, NullSpace,NODES,ELEMENTS] = BuildElasticity(NELS, varargin)

% Get the matrix
[MAT,NODES,ELEMENTS]=elasticfun(NELS,varargin{:});

% Get the RBM
NullSpace=build_elastic_rbm(NODES);

% Build Amat
newMap=Map(size(MAT,1),1);
Amat = Operator(MAT,newMap,newMap,@MatlabApply);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RBM=build_elastic_rbm(NODES)
% Constants
[N,dim]=size(NODES);
NDOF=[1,3,6];

% Allocs
RBM=zeros(dim*N,NDOF(dim));

% Translational DOFs / INDICES
for I=1:dim,
  IDX{I}=I:dim:dim*N;
  RBM(IDX{I},I)=ones(N,1);
end

% Recenter nodes
CTR=sum(NODES) / N;
CNODES = NODES - repmat(CTR,N,1);

% Rotational DOF:  Thanks to Farhat, Pierson and Lesoinne 2000,
% equation (52), you know, once I've managed to make sense out of
% the blasted thing.
if(dim>=2)
  % Rotate in X-Y Plane (around Z axis): [-y ;x];
  RBM(IDX{1},dim+1)=-CNODES(:,2);
  RBM(IDX{2},dim+1)= CNODES(:,1);
end
if(dim==3)
  % Rotate in Y-Z Plane (around X axis): [-z;y]
  RBM(IDX{2},dim+2)=-CNODES(:,3);
  RBM(IDX{3},dim+2)= CNODES(:,2);

  % Rotate in X-Z Plane (around Y axis): [z ;-x]
  RBM(IDX{1},dim+3)= CNODES(:,3);
  RBM(IDX{3},dim+3)=-CNODES(:,1);
end

RBM=normalize(RBM);


% Happy normalize function (normalizes columns)
function [AN,NF]=normalize(A)
SZ=size(A,2);
NF=zeros(SZ,1);
AN=0*A;
for I=1:SZ,
  NF(I)=norm(A(:,I));
  if(NF(I) < 1e-11)
    NF(I)=1;AN(:,I)=A(:,I);
  else
   AN(:,I)=A(:,I)/NF(I);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [MAT,NODES]=elasticfun(NELS,[E],[nu],[stretch],[mode])
%
% Generates a discretization of linear elasticity in two dimensions
% under either plane stress or plain strain assumptions or in three
% dimensions.  Dirichlet conditions are i
% Input:
% NELS    - Number of elements in each dimension
% E       - Elastic modulus (or function to compute it)
% nu      - Poisson's ration
% stretch - Mesh stretching
% mode    - [2D only] - 'plane stress' or 'plane strain'
%
%
% Output:
% MAT     - Discretizes operator
% NODES   - Location of the nodes
%
% by: Chris Siefert <csiefer@sandia.gov>
% Last Update: 05/05/10 <csiefer>
%
function [MAT,NODES,ELEMENTS]=elasticfun(NELS,varargin)
% Check for E, nu
if(nargin>1), E=varargin{1}; else E=1e9;end
if(nargin>2), nu=varargin{2}; else nu=.25;end
if(nargin>3), stretch=varargin{3}; else stretch=[1,1,1];end
if(nargin>4), mode=varargin{4}; else mode='plane stress';end

% Constants
dim=length(NELS);

% Mesh the problem
[NODES,ELEMENTS,BIDX,BNODES]=mesher(NELS,stretch);

% Build the stiffness matrix
if(dim==2), MAT=elasticfun_2d(NODES,ELEMENTS,E,nu,mode);
elseif(dim==3), MAT=elasticfun_3d(NODES,ELEMENTS,E,nu,mode); end

% Apply boundary conditions
IIDX=setdiff(1:size(MAT,1),BIDX);
MAT=MAT(IIDX,IIDX);
NODE_REINDEX=0*NODES;
NIDX=setdiff(1:size(NODES,1),BNODES);
NODES=NODES(NIDX,:);
NODE_REINDEX(NIDX)=1:length(NIDX);

% Blast all boundary-containing elements
NEW_ELS=[];
for I=1:size(ELEMENTS,1),
  ridx=NODE_REINDEX(ELEMENTS(I,:));
  ii=find(ridx==0);
  if(length(ii)==0),
    NEW_ELS=[NEW_ELS;ridx];
  end
end
ELEMENTS=NEW_ELS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NODES,ELEMENTS,BIDX,BDY]=mesher(NELS,stretch)
% Constants
NPTS=NELS+1;
Nn=prod(NPTS);
Nel=prod(NELS);
dim=length(NELS);

ROOT_NODES{1}=[1:NPTS(1)]'*stretch(1);
ROOT_NODES{2}=[1:NPTS(2)]'*stretch(2);
if(dim==3), ROOT_NODES{3}=[1:NPTS(3)]'*stretch(3);end

% Nodal Locations (2D + 3D)
NODES=ROOT_NODES{1};
for I=2:dim,
  SZ=size(NODES,1);
  NODES=[repmat(NODES,NPTS(I),1),reshape(repmat(ROOT_NODES{I}',SZ,1),SZ*NPTS(I),1)];
end

% Get the nodes in the bottom left corner of the elements
IDX=coord_to_idx(idx_to_coord([1:Nel]',NELS),NPTS);

if(dim==2),
  % Element matrix (2D) Ordering
  %  4-3
  %  | |
  %  1-2
  %
  ELEMENTS=[IDX,IDX+1,IDX+NPTS(1)+1,IDX+NPTS(1)];

  % Find ID's of boundary nodes
  %BDY=find(NODES(:,1)==stretch(1) | NODES(:,1)==NPTS(1)*stretch(1) | ...
  %           NODES(:,2)==stretch(2) | NODES(:,2)==NPTS(2)*stretch(2));

  BDY=find(NODES(:,1)==stretch(1));
elseif(dim==3),
  % Element matrix (3D) Ordering
  %     4-3           8-7
  %     1-2           5-6
  %  Bottom Plane   Top Plane
  % (see exodus node ordering)
  ELEMENTS(:,1:4)=[IDX,IDX+1,IDX+NPTS(1)+1,IDX+NPTS(1)];
  ELEMENTS(:,5:8)=ELEMENTS(:,1:4)+NPTS(1)*NPTS(2);

  % Find ID's of boundary nodes
  %  BDY=find(NODES(:,1)==stretch(1) | NODES(:,1)==NPTS(1)*stretch(1) | ...
  %         NODES(:,2)==stretch(2) | NODES(:,2)==NPTS(2)*stretch(2) | ...
  %         NODES(:,3)==stretch(3) | NODES(:,3)==NPTS(3)*stretch(3));

  BDY=find(NODES(:,1)==stretch(1));
end

% Nodes to DOF
BIDX=reshape(repmat((BDY-1)*dim,1,dim)'+repmat(1:dim,length(BDY),1)',dim*length(BDY),1);

%BIDX=[];
%BDY=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function COORD=idx_to_coord(IDX,NITM)
dim=length(NITM);
COORD(:,1)=1+mod(IDX-1,NITM(1));
if(dim==2),COORD(:,2)=1+(IDX-COORD(:,1)) / NITM(1);end
if(dim==3),
  COORD(:,2)=1+mod((IDX-COORD(:,1))/NITM(1),NITM(2));
  COORD(:,3)=1+(IDX-COORD(:,1)-(COORD(:,2)-1)*NITM(1))/(NITM(1)*NITM(2));
end

function IDX=coord_to_idx(COORD,NITM)
dim=length(NITM);
if(dim==2), IDX=1+(COORD(:,1)-1) + (COORD(:,2)-1)*NITM(1); end
if(dim==3), IDX=1+(COORD(:,1)-1) + (COORD(:,2)-1)*NITM(1) + (COORD(:,3)-1)*NITM(1)*NITM(2); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K=elasticfun_2d(NODES,ELEMENTS,E,nu,mode)
% Constants
Nel=size(ELEMENTS,1);
Nn=size(NODES,1);
dim=2;

% Reference Element Coords
rcoord=[-1,-1; 1,-1; 1,1; -1,1];
NperEl=4;

% Condensed DOF
CDOF=8;

% Material constants
t=1;

% Build Material Matrix
if(strcmp(mode,'plane stress'))
  if(isnumeric(E)), D=E / (1-nu^2) * [1,nu,0;nu,1,0;0,0,(1-nu)/2];
  else D0=1 / (1-nu^2) * [1,nu,0;nu,1,0;0,0,(1-nu)/2]; end
elseif(strcmp(mode,'plane strain'))
  if(isnumeric(E)), D=E / (1+nu)/(1-2*nu) * [1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
  else D0=1 / (1+nu)/(1-2*nu) * [1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2]; end
else
  fprintf('Error: Unknown Material Model for 2D\n');
  return;
end


% Gauss Points (reference)
ORT=1/sqrt(3);
GP=[ORT,ORT;ORT,-ORT;-ORT,ORT;-ORT,-ORT];
NGP=size(GP,1);
WT=ones(1,NGP);

% Evaluate the B matrix for the reference element
B=zeros(3,CDOF);
for J=1:NGP,
  dxi=eval2_Ni_dxi(rcoord,GP(J,1),GP(J,2));
  deta=eval2_Ni_deta(rcoord,GP(J,1),GP(J,2));
  B(1,1:2:7) = dxi(1:4);
  B(2,1:2:7) = deta(1:4);
  B(3,2:2:8) = dxi(1:4);
  B(4,2:2:8) = deta(1:4);
  BGP{J}=B;
  SFD{J}=[dxi;deta];
end

% See 6.2-9 from Cook
REORDER=[1,0,0,0;...
         0,0,0,1;...
         0,1,1,0];


% Evaluate each gauss point the stiffness matrix
K=sparse(dim*Nn,dim*Nn);
for I=1:Nel,
  NIDX=ELEMENTS(I,:);
  DIDX=node_to_dof(NIDX',dim);
  KE=sparse(CDOF,CDOF);

  if(~isnumeric(E)), D=D0*feval(E,mean(NODES(ELEMENTS(I,:),:)));end

  for J=1:NGP,
    B=BGP{J};
    SF=SFD{J};
    JAC=SF(:,1:NperEl)*NODES(NIDX,:);
    detJ=det(JAC);
    J2=[JAC,zeros(2,2);zeros(2,2),JAC];
    J2B=REORDER*(J2\BGP{J});
    KE = KE  + WT(J) *t*J2B'*D *J2B* detJ;
  end
  K=ssubsasgn(K,DIDX,DIDX,KE);
end

% Hard Symmetry, cleanup zeros
K=(K+K')/2;
K=K .* (abs(K)>=1e-15*norm(K,inf));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv=eval2_Ni_dxi(rcoord, xi, eta)
rv(1:4)=[rcoord(:,1).*(1 + eta*rcoord(:,2))/4]';
rv(5)=-2*xi;
rv(6)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv=eval2_Ni_deta(rcoord, xi, eta)
rv(1:4)=[(1 + xi*rcoord(:,1)).*rcoord(:,2)/4]';
rv(5)=0;
rv(6)=-2*eta;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=elasticfun_3d(NODES,ELEMENTS,E,nu,mode)
% Constants
Nel=size(ELEMENTS,1);
Nn=size(NODES,1);
dim=3;

% Reference Element Coords
rcoord=[ -1,-1,-1;
         1,-1,-1;
         1, 1,-1;
        -1, 1,-1;
        -1,-1, 1;
         1,-1, 1;
         1, 1, 1;
        -1, 1, 1];
NperEl=8;

% Uncondensed DOF
%
% Condensed DOF
%
DOF=24;

% Build Material Matrix
if(isnumeric(E)), c=E / ((1+nu)*(1-2*nu));
else c=1 / ((1+nu)*(1-2*nu)); end
M1=[1-nu,nu,nu;...
    nu,1-nu,nu;...
    nu,nu,1-nu];
D=c*[M1,zeros(3);zeros(3),(1-2*nu)/2*eye(3)];
D0=D;

% Materials
t=1;

% Gauss Points (reference)
ORT=1/sqrt(3);
GP=[ ORT, ORT, ORT;...
     ORT,-ORT, ORT;...
    -ORT, ORT, ORT;...
    -ORT,-ORT, ORT;...
     ORT, ORT,-ORT;...
     ORT,-ORT,-ORT;...
    -ORT, ORT,-ORT;...
    -ORT,-ORT,-ORT];
NGP=size(GP,1);
WT=ones(1,NGP);

% Evaluate the B matrix for the reference element
for J=1:NGP,
  B2=zeros(9,DOF);
  dxi  =eval3_Ni_dxi(rcoord,GP(J,1),GP(J,2),GP(J,3));
  deta =eval3_Ni_deta(rcoord,GP(J,1),GP(J,2),GP(J,3));
  dzeta=eval3_Ni_dzeta(rcoord,GP(J,1),GP(J,2),GP(J,3));

  B2(1,1:3:24)=dxi;
  B2(2,1:3:24)=deta;
  B2(3,1:3:24)=dzeta;
  B2(4,2:3:24)=dxi;
  B2(5,2:3:24)=deta;
  B2(6,2:3:24)=dzeta;
  B2(7,3:3:24)=dxi;
  B2(8,3:3:24)=deta;
  B2(9,3:3:24)=dzeta;

  BGP2{J}=B2;
  SFD{J}=[dxi;deta;dzeta];
end

% See 6.5-3 from Cook
REORDER=[1,0,0,0,0,0,0,0,0;...
         0,0,0,0,1,0,0,0,0;...
         0,0,0,0,0,0,0,0,1;...
         0,1,0,1,0,0,0,0,0;...
         0,0,0,0,0,1,0,1,0;...
         0,0,1,0,0,0,1,0,0];

% Evaluate each gauss point the stiffness matrix
K=sparse(dim*Nn,dim*Nn);
for I=1:Nel,
  NIDX=ELEMENTS(I,:);
  DIDX=node_to_dof(NIDX',dim);
  KE=sparse(DOF,DOF);
  if(~isnumeric(E)), D=D0*feval(E,mean(NODES(ELEMENTS(I,:),:)));end

  for J=1:NGP,
    SF=SFD{J};
    JAC=SF*NODES(NIDX,:);
    detJ=det(JAC);
    J2=[JAC,zeros(3,6);zeros(3,3),JAC,zeros(3,3);zeros(3,6),JAC];
    J2B=REORDER*(J2\BGP2{J});
    KE = KE  + WT(J) *t*J2B'*D *J2B* detJ;
  end
  K(DIDX,DIDX) = K(DIDX,DIDX) + KE;
  %  K=ssubsasgn(K,DIDX,DIDX,KE(1:DOF,1:DOF));
end

% Hard Symmetry, cleanup zeros
K=(K+K')/2;
K=K .* (abs(K)>=1e-15*norm(K,inf));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv=eval3_Ni_dxi(rcoord,xi,eta,zeta)
rv=[rcoord(:,1).*(1+rcoord(:,2)*eta).*(1+rcoord(:,3)*zeta) / 8]';
function rv=eval3_Ni_deta(rcoord,xi,eta,zeta)
rv=[rcoord(:,2).*(1+rcoord(:,1)*xi ).*(1+rcoord(:,3)*zeta) / 8]';
function rv=eval3_Ni_dzeta(rcoord,xi,eta,zeta)
rv=[rcoord(:,3).*(1+rcoord(:,1)*xi ).*(1+rcoord(:,2)*eta ) / 8]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function DIDX=node_to_dof(IDX,ndof)
% Converts from nodal indices to DOF indices, assuming that there
% are ndof DOFs per node.  Assumes that numbering is done by node,
% with the DOFs numbered in order.
% Input:
% IDX     - Nodal indices
% ndof    - Number of DOFs per node
% Output:
% DIDX    - DOF indices corresponding to nodes IDX
%
% By: Chris Siefert <csiefer@sandia.gov>
% Last Updated: 08/22/2006 <csiefer>

function DIDX=node_to_dof(IDX,ndof)
DIDX=reshape(repmat(ndof*(IDX'-1)+1,ndof,1) + repmat([0:ndof-1]',1,length(IDX)),ndof*length(IDX),1);

