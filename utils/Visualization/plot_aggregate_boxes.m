function h=plot_aggregate_boxes(NODES,agg,dx,varargin)
if(nargin==4), BCOLOR=varargin{1}; else BCOLOR='b';end

Na=max(agg);
dim=size(NODES,2);
if(dim~=2 && dim~=3), return; end
h=[];

for I=1:Na,
  IDX=find(agg==I);
  AMAX=max(NODES(IDX,:)) + dx;
  AMIN=min(NODES(IDX,:)) - dx;
  
  if(dim==2), 
    h=[h,plot_hull_2d(NODES(IDX,:),BCOLOR)]; 
  elseif(dim==3), h=[h,plot_box_3d(AMIN,AMAX,BCOLOR)]; end


end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots a convex hull in 2D
function h=plot_hull_2d(RNODES,BCOLOR)
Ni=size(RNODES,1);
CENTER=mean(RNODES);
VN=1.2*RNODES-.2*repmat(CENTER,Ni,1);
K=convhulln(VN);

X=[VN(K(:,1),1);VN(K(1,1),1)];
Y=[VN(K(:,1),2);VN(K(1,1),2)];
h=patch(X,Y,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots a box in 2D (deprecated)
function h=plot_box_2d(C1,C2,BCOLOR)
X=[C1(:,1),C2(:,1),C2(:,1),C1(:,1),C1(:,1)];
Y=[C1(:,2),C1(:,2),C2(:,2),C2(:,2),C1(:,2)];
h=patch(X,Y,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots a box in 3D
function h=plot_box_3d(C1,C2,BCOLOR)
h=[];

% X-Y planes
X=[C1(:,1),C2(:,1),C2(:,1),C1(:,1),C1(:,1)];
Y=[C1(:,2),C1(:,2),C2(:,2),C2(:,2),C1(:,2)];
Z=repmat(C1(:,3),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];
Z=repmat(C2(:,3),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];

% X-Z planes
X=[C1(:,1),C2(:,1),C2(:,1),C1(:,1),C1(:,1)];
Z=[C1(:,3),C1(:,3),C2(:,3),C2(:,3),C1(:,3)];
Y=repmat(C1(:,2),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];
Y=repmat(C2(:,2),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];

% Y-Z planes
Y=[C1(:,2),C2(:,2),C2(:,2),C1(:,2),C1(:,2)];
Z=[C1(:,3),C1(:,3),C2(:,3),C2(:,3),C1(:,3)];
X=repmat(C1(:,1),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];
X=repmat(C2(:,1),5,1);
h=[h,patch(X,Y,Z,BCOLOR,'FaceAlpha',0,'FaceColor',[1,1,1],'EdgeColor',BCOLOR)];


