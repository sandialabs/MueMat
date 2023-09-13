% Copyright 2010  Institute for Computational Mechanics, TUM
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library.  If not, see doc/GNU_Lesser_GPL.txt.

classdef Scatra2D_vg1
   properties (Access = private)
       problem_ = [];

       gausspt_ = [];
       shpfct_ = [];
       dis_ = [];
   end

   methods (Access = public)
       function [this] = Scatra2D_vg1(nele1,nele2,kappa,frhs,afun,description)

           %load geometry.mat      % cube [0,1]^2
           %[p,e,t] = Scatra2D.poi_mesh(g,nele1,nele2);
           load mesh32.mat

           % gauss points
           this.gausspt_.nquad = 3;
           this.gausspt_.wgt = [1/6; 1/6; 1/6];
           this.gausspt_.xg = [1/6 1/6; 2/3 1/6; 1/6 2/3];

           % shape functions
           this.shpfct_.f1 = inline('1-r-s','r','s');
           this.shpfct_.f2 = inline('r','r','s');
           this.shpfct_.f3 = inline('s','r','s');
           this.shpfct_.f1d = inline('[-1;-1]','r','s');
           this.shpfct_.f2d = inline('[1;0]','r','s');
           this.shpfct_.f3d = inline('[0;1]','r','s');

           this.dis_.nElements = size(t,2);

           for iele=1:this.dis_.nElements
               temp = p(:,t(1:3,iele));
               ele.xyze = zeros(3,2);
               ele.xyze(1,:) = temp(:,1)';  % real coords (node 1)
               ele.xyze(2,:) = temp(:,2)';  % real coords (node 2)
               ele.xyze(3,:) = temp(:,3)';  % real coords (node 3)
               this.dis_.e(iele) = ele;
               clear temp;
               clear ele;
           end
           clear iele;

           numdofpernode = 1;
           numnodes = size(p(1,:),2);
           numdofs = numnodes * numdofpernode;

           this.problem_.dis = this.dis_;
           this.problem_.shpfct = this.shpfct_;
           this.problem_.gausspt = this.gausspt_;
           this.problem_.geo.p = p;
           this.problem_.geo.e = e;
           this.problem_.geo.t = t;
           this.problem_.numnodes = numnodes;
           this.problem_.numdofpernode = numdofpernode;
           this.problem_.numdofs = numdofs;

           % TODO: make bcs accesible
%            this.problem_.bc.dirichnodes = find(p(1,:)==0 | p(2,:)==0 | p(1,:)==1);             % left, bottom and right border
%            this.problem_.bc.dirichvals = zeros(length(this.problem_.bc.dirichnodes)*numdofpernode,1);
%            this.problem_.bc.dirichnodes = [this.problem_.bc.dirichnodes, find(p(2,:)==1)];             % top border
%            this.problem_.bc.dirichvals = [this.problem_.bc.dirichvals;ones(length(find(p(2,:)==1)),1)];
           this.problem_.bc.dirichnodes = find((p(1,:)==0 & p(2,:) < 0.7) | p(2,:)==0 | p(1,:)==1);             % left, bottom and right border
           this.problem_.bc.dirichvals = zeros(length(this.problem_.bc.dirichnodes)*numdofpernode,1);
           this.problem_.bc.dirichnodes = [this.problem_.bc.dirichnodes, find(p(2,:)==1 | (p(1,:)==0)&p(2,:)>=0.7)];             % top border
           this.problem_.bc.dirichvals = [this.problem_.bc.dirichvals;ones(length(find(p(2,:)==1 | (p(1,:)==0)&p(2,:)>=0.7)),1)];

           if varexist('frhs') this.problem_.frhs = frhs;
           else this.problem_.frhs = inline('[0;0;0]','x','y'); end;

           if varexist('afun') this.problem_.afun = afun;
           %else this.problem_.afun = inline('[cos(pi/6),sin(pi/6)]','x','y'); end;
           else this.problem_.afun = inline('[cos(-pi/3),sin(-pi/3)]','x','y'); end;

           if varexist('kappa') this.problem_.kappa = kappa;
           %else this.problem_.kappa = 1/3200 * ones(this.problem_.numdofpernode,1); end;
           else this.problem_.kappa = 1e-6 * ones(this.problem_.numdofpernode,1); end;

           if varexist('description') this.problem_.description = description;
           else this.problem_.description = 'Scatra 2d doubleglazing'; end;

           clear numdofpernode; clear numnodes; clear numdofs; clear h; clear xnod; clear i; clear nele;
           clear p; clear e; clear t;
       end

       function [A,b,problem] = Build(this)
           K = this.assemble_K();
           C = this.assemble_C();
           F = this.assemble_F();

           A = K + C;
           b = F;

           newMap = Map(size(A,1), 1); % generate full map (point-wise)
           A = Operator(A,newMap,newMap,@MatlabApply, this.problem_.description);
           clear newMap;

           problem = this.problem_;

       end

       function plot_vector(this,vector)
           Scatra2D.plot_scatra2d(this.problem_.geo.p,this.problem_.geo.e,this.problem_.geo.t,'zdata',vector);
       end

%        function visualizeAggsFromP(this,PP)
%            hold on
%            coords = this.problem_.geo.p;
%
%            % loop over all aggregates (= cols of PP)
%            for j=1:size(PP,2)
%                fcol = [0.5,0.5+0.5*rand(1,1),0.5*rand(1,1)];
%                nodeIds = find(PP(:,j));
%                aggCoords = coords(:,nodeIds);
%
%                if(size(aggCoords,2) > 1)
%                    K = convhull(aggCoords(1,:),aggCoords(2,:));
%                    alpha(0.1)
%                    fill(aggCoords(1,K),aggCoords(2,K),fcol);
%                else
%                    plot(aggCoords(1,:),aggCoords(2,:),'rx');
%                    alpha(0.5);
%                end
%            end
%            hold off
%        end

   end

   methods (Static = true)
       [p,r,t] = poi_mesh(g,n1,n2);
       [x,y]=pde_igeom(dl,bs,s);
       h=plot_scatra2d(p,e,t,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14);
       [uxy,tn,al2,al3]=tri2grid(p,t,u,tn,al2,al3);
   end

   methods (Access = private)
       function [K] = assemble_K(this)

           K = sparse(this.problem_.numdofs,this.problem_.numdofs);

           nel = this.problem_.dis.nElements;


           for ele=1:nel

               iel = 3;    % Knoten pro Element
               kappa = this.problem_.kappa; %ones(this.problem_.numdofpernode,1);
               estif = zeros(iel*this.problem_.numdofpernode);

               % loop over all gauss pts
               for iquad=1:this.problem_.gausspt.nquad

                   [funct,derxy,fac] = this.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,this.problem_.gausspt,this.problem_.shpfct,this.problem_.dis);


                   for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                       for ui=1:iel   % alle Knoten
                           for vi=1:iel  % alle Testfunktionen
                               estif((ui-1)*this.problem_.numdofpernode+dofindex,(vi-1)*this.problem_.numdofpernode+dofindex) = ...
                                   estif((ui-1)*this.problem_.numdofpernode+dofindex,(vi-1)*this.problem_.numdofpernode+dofindex) + ...
                                   fac * kappa(dofindex) * ...
                                   (derxy(ui,1)*derxy(vi,1) + derxy(ui,2)*derxy(vi,2));
                               % N_x*M_x + N_y*M_y = (\nabla N, \nabla M)
                               % d.h. auf T_0 normierte Variante von (\nabla Phi,
                               % \nabla Psi)
                           end
                       end
                   end
                   %estif % Elementsteifigkeitsmatrix f�r ele=1 und Gau�punkt iquad
               end

               for dofindex=1:this.problem_.numdofpernode
                   % node numbers of element
                   nodes = this.problem_.geo.t(1:3,ele);
                   K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,dofindex);
                   K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,1*this.problem_.numdofpernode+dofindex);
                   K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,2*this.problem_.numdofpernode+dofindex);
                   K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,dofindex);
                   K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,1*this.problem_.numdofpernode+dofindex);
                   K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,2*this.problem_.numdofpernode+dofindex);
                   K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,dofindex);
                   K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,1*this.problem_.numdofpernode+dofindex);
                   K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = K((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,2*this.problem_.numdofpernode+dofindex);
               end
           end

           % Dirichlet RB
           bcnodes = this.problem_.bc.dirichnodes;

           for bc = 1 : size(bcnodes,2)
               for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                   K((bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex,:) = 0;
                   K((bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex,(bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex) = 0.5;
               end
           end
       end % end assemble_K

       function [C] = assemble_C(this)

           C = sparse(this.problem_.numdofs,this.problem_.numdofs);

           nel = this.problem_.dis.nElements;


           for ele=1:nel

               iel = 3;    % nodes per element
               kappa = this.problem_.kappa;
               estif = zeros(iel*this.problem_.numdofpernode);

               % nodes of current element
               nodes = this.problem_.geo.t(1:3,ele);

               % calculate stabilization parameter
               tau = this.calc_tau(this.problem_,ele);

               % loop over all gauss points
               for iquad=1:this.problem_.gausspt.nquad

                   % shape functions and shape derivatives
                   % at current gauss point
                   [funct,derxy,fac] = this.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,this.problem_.gausspt,this.problem_.shpfct,this.problem_.dis);

                   % node coordinates
                   coords1 = this.problem_.geo.p(:,nodes(1));
                   coords2 = this.problem_.geo.p(:,nodes(2));
                   coords3 = this.problem_.geo.p(:,nodes(3));
                   xgausp = zeros(this.problem_.numdofpernode,1);
                   xgausp = coords1 + this.problem_.gausspt.xg(iquad,1) * (coords2-coords1) + this.problem_.gausspt.xg(iquad,2) * (coords3-coords1);


                   for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                       for vi=1:iel  % alle Testfunktionen
                           for ui=1:iel   % alle dofs
                               %                         % Koordinaten
                               %                                                 nodes = this.problem_.geo.t(1:3,ele);
                               %                                                 coord = this.problem_.geo.p(:,nodes(ui));
                               %                                                 a = zeros(1,2);
                               %                                                 a = this.problem_.afun(coord(1),coord(2));
                               %                         %                         this.problem_.a = a;
                               a = zeros(1,2);
                               a = this.problem_.afun(xgausp(1),xgausp(2));

                               % (a, \nabla u) * v
                               estif((vi-1)*this.problem_.numdofpernode+dofindex,(ui-1)*this.problem_.numdofpernode+dofindex) = ...
                                   estif((vi-1)*this.problem_.numdofpernode+dofindex,(ui-1)*this.problem_.numdofpernode+dofindex) + ...
                                   fac * (a(1)*derxy(ui,1) + a(2)*derxy(ui,2)) * funct(vi);

                               % stabilization term
                               % only for stationary this.problem_s!
                               estif((vi-1)*this.problem_.numdofpernode+dofindex,(ui-1)*this.problem_.numdofpernode+dofindex) = ...
                                   estif((vi-1)*this.problem_.numdofpernode+dofindex,(ui-1)*this.problem_.numdofpernode+dofindex) + ...
                                   tau(dofindex) * fac * (a(1)*derxy(ui,1) + a(2)*derxy(ui,2)) * (a(1)*derxy(vi,1) + a(2)*derxy(vi,2));
                           end
                       end
                   end
               end

               for dofindex=1:this.problem_.numdofpernode
                   % node numbers of element
                   nodes = this.problem_.geo.t(1:3,ele);
                   C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,dofindex);
                   C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,1*this.problem_.numdofpernode+dofindex);
                   C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(1)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(dofindex,2*this.problem_.numdofpernode+dofindex);
                   C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,dofindex);
                   C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,1*this.problem_.numdofpernode+dofindex);
                   C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(2)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(1*this.problem_.numdofpernode+dofindex,2*this.problem_.numdofpernode+dofindex);
                   C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(1)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,dofindex);
                   C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(2)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,1*this.problem_.numdofpernode+dofindex);
                   C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) = C((nodes(3)-1)*this.problem_.numdofpernode+dofindex,(nodes(3)-1)*this.problem_.numdofpernode+dofindex) + estif(2*this.problem_.numdofpernode+dofindex,2*this.problem_.numdofpernode+dofindex);
               end
           end

           % Dirichlet RB
           bcnodes = this.problem_.bc.dirichnodes;

           for bc = 1 : size(bcnodes,2)
               for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                   C((bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex,:) = 0;
                   C((bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex,(bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex) = 0.5;   % here we should have 1, but we assemble a 1 within setup routine for K
               end
           end
       end % end assemble_C

       function [F] = assemble_F(this)

           iel = 3;    % nodes per element
           nel = this.problem_.dis.nElements;

           % rhs
           F = zeros(this.problem_.numdofs,1);
           for ele = 1:nel
               % nodes of current element
               nodes = this.problem_.geo.t(1:3,ele);

               f = zeros(iel*this.problem_.numdofpernode,1);

               % Stabilisierungsparameter bestimmen
               tau = this.calc_tau(this.problem_,ele);


               % loop over all gaus pts
               for iquad = 1:this.problem_.gausspt.nquad

                   [funct,derxy,fac] = this.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,this.problem_.gausspt,this.problem_.shpfct,this.problem_.dis);

                   % node coords
                   coords1 = this.problem_.geo.p(:,nodes(1));
                   coords2 = this.problem_.geo.p(:,nodes(2));
                   coords3 = this.problem_.geo.p(:,nodes(3));
                   xgausp = zeros(this.problem_.numdofpernode,1);
                   xgausp = coords1 + this.problem_.gausspt.xg(iquad,1) * (coords2-coords1) + this.problem_.gausspt.xg(iquad,2) * (coords3-coords1);

                   r = zeros(iel*this.problem_.numdofpernode,1);
                   r = this.problem_.frhs(xgausp(1),xgausp(2));


                   for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                       for vi=1:iel  % alle Testfunktionen
                           f((vi-1)*this.problem_.numdofpernode+dofindex,1) = f((vi-1)*this.problem_.numdofpernode+dofindex,1) + fac*funct(vi) * r((vi-1)*this.problem_.numdofpernode+dofindex,1);

                           %                     % Konvektionsgeschwindigkeit für aktuellen Punkt
                           %                     coord = this.problem_.geo.p(:,nodes(vi));
                           %                     a = zeros(1,2);
                           %                     a = this.problem_.afun(coord(1),coord(2));
                           a = zeros(1,2);
                           a = this.problem_.afun(xgausp(1),xgausp(2));

                           % Stabilisierungsparameter
                           % konvektiv taufac * conv(vi) * rhsint
                           f((vi-1)*this.problem_.numdofpernode+dofindex,1) = f((vi-1)*this.problem_.numdofpernode+dofindex,1) + ...
                               tau(dofindex) * fac * (a(1)*derxy(vi,1) + a(2)*derxy(vi,2)) * r((vi-1)*this.problem_.numdofpernode+dofindex,1);
                       end
                       %f
                   end
               end

               % assembly
               for dofindex=1:this.problem_.numdofpernode
                   F((nodes(1)-1)*this.problem_.numdofpernode+dofindex,1) = F((nodes(1)-1)*this.problem_.numdofpernode+dofindex,1) + f(dofindex,1);
                   F((nodes(2)-1)*this.problem_.numdofpernode+dofindex,1) = F((nodes(2)-1)*this.problem_.numdofpernode+dofindex,1) + f(1*this.problem_.numdofpernode+dofindex,1);
                   F((nodes(3)-1)*this.problem_.numdofpernode+dofindex,1) = F((nodes(3)-1)*this.problem_.numdofpernode+dofindex,1) + f(1*this.problem_.numdofpernode+dofindex,1);
               end
           end

           % Dirichlet RB
           bcnodes = this.problem_.bc.dirichnodes;

           for bc = 1 : size(bcnodes,2)
               for dofindex=1:this.problem_.numdofpernode % alle dofindizes
                   F((bcnodes(bc)-1)*this.problem_.numdofpernode+dofindex,1)=this.problem_.bc.dirichvals((bc-1)*this.problem_.numdofpernode+dofindex);
               end
           end

       end % end assemble_F

       % current Gausspt: iquad (1-3)
       % current element: ele
       % derxy: global derivatives (3x2 matrix)
       function [funct,derxy,fac] = EvalShapeFuncAndDerivAtIntPoint(this,iquad,ele,gausspt,shpfct,dis)
           % coords of current gauss pt
           e1 = gausspt.xg(iquad,1);
           e2 = gausspt.xg(iquad,2);

           % evaluate the three shape functions at gauss pt (e1,e2)
           funct = zeros(3,1);
           funct(1) = shpfct.f1(e1,e2);
           funct(2) = shpfct.f2(e1,e2);
           funct(3) = shpfct.f3(e1,e2);

           % calculate the three derivatives of the shape functions at gauss pt (e1,e2)
           deriv = zeros(3,2);
           deriv(1,:) = shpfct.f1d(e1,e2)';
           deriv(2,:) = shpfct.f2d(e1,e2)';
           deriv(3,:) = shpfct.f3d(e1,e2)';

           iel = 3;    % number of points per element

           xjm = zeros(2,2);

           %          +-            -+ T      +-            -+
           %          | dx   dx   dx |        | dx   dy   dz |
           %          | --   --   -- |        | --   --   -- |
           %          | dr   ds   dt |        | dr   dr   dr |
           %          |              |        |              |
           %          | dy   dy   dy |        | dx   dy   dz |
           %  xjm  =  | --   --   -- |   =    | --   --   -- |
           %          | dr   ds   dt |        | ds   ds   ds |
           %          |              |        |              |
           %          | dz   dz   dz |        | dx   dy   dz |
           %          | --   --   -- |        | --   --   -- |
           %          | dr   ds   dt |        | dt   dt   dt |
           %          +-            -+        +-            -+
           % here of course only 2d!
           for i=1:2
               for j=1:2
                   dum = 0;
                   for l=1:iel
                       dum = dum + deriv(l,i)*dis.e(ele).xyze(l,j);
                   end
                   xjm(i,j) = dum;
               end
           end

           % determinant of xjm due to Sarrus' formula (~> det (J) )
           determinante = xjm(1,1)*xjm(2,2) - xjm(1,2)*xjm(2,1);

           % check det < 0?
           if determinante < 0
               error('det < 0');
           end

           % gauss weights * det(J)
           fac = gausspt.wgt(iquad)*determinante;

           % inverse of xjm
           xij = zeros(2,2);
           xij(1,1) = +xjm(2,2)/determinante;
           xij(2,1) = -xjm(2,1)/determinante;
           xij(1,2) = -xjm(1,2)/determinante;
           xij(2,2) = +xjm(1,1)/determinante;

           % global derivatives
           derxy = zeros(3,2);
           for k=1:iel
               derxy(k,1) = xij(1,1)*deriv(k,1) + xij(1,2)*deriv(k,2);
               derxy(k,2) = xij(2,1)*deriv(k,1) + xij(2,2)*deriv(k,2);
           end

       end % EvalShapeFuncAndDerivAtIntPoint

       % stabilization for scatra due to Franca and Valentin 2000
       function [tau] = calc_tau(this,problem,cur_element)

           % evaluate shape functions
           e1 = 1/3;
           e2 = 1/3;
           funct = zeros(3,1);
           funct(1) = problem.shpfct.f1(e1,e2);
           funct(2) = problem.shpfct.f2(e1,e2);
           funct(3) = problem.shpfct.f3(e1,e2);

           % area of triangle with Heron's formula
           xyze = problem.dis.e(cur_element).xyze;  % coords of current element (3x2)
           a = (xyze(1,1)-xyze(2,1))^2 + (xyze(1,2)-xyze(2,2))^2;
           b = (xyze(2,1)-xyze(3,1))^2 + (xyze(2,2)-xyze(3,2))^2;
           c = (xyze(3,1)-xyze(1,1))^2 + (xyze(3,2)-xyze(1,2))^2;
           area = 0.25 * sqrt(2*a*b+2*b*c+2*c*a-a^2-b^2-c^2);

           hk = sqrt(area);        % h-const for current element

           mk = 1/3;               % tau constant

           velint = zeros(2,1);    % 2 dim velocity vector at element center

           iel = 3;        % nodes per element
           for j=1:iel
               nodes = problem.geo.t(1:3,cur_element);
               coord = problem.geo.p(:,nodes(j));
               a = zeros(1,2);
               a = problem.afun(coord(1),coord(2));
               velint(1,1) = velint(1,1) + funct(j) * a(1);    % hier eigentlich evel(i+(2*j)) mit i=1
               velint(2,1) = velint(2,1) + funct(j) * a(2);    % hier eigentlich evel(i+2(*j)) mit i=2
           end

           vel_norm = norm(velint);

           % stabilization parameter for stationary case
           tau = zeros(problem.numdofpernode,1);
           for k=1:problem.numdofpernode
               epe2 = mk * vel_norm * hk / problem.kappa(k);
               xi2 = max(epe2,1.0);
               tau(k) = ((hk^2) * mk)/(2*problem.kappa(k)*xi2);
           end

       end

   end% end methods (private)

   methods (Access = protected)
       function Copy_(this,src,mc)
           [cmd,data,mc] = this.CopyCmd_(src,mc);
           eval(cmd);
       end
   end % protected methods
end