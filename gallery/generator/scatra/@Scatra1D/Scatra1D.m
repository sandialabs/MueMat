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

classdef Scatra1D
   properties (Access = private)
       problem_ = [];

       xstart_ = 0;
       xend_ = 1;

       gausspt_ = [];
       shpfct_ = [];
       dis_ = [];
   end

   methods (Access = public)
       function [this] = Scatra1D(nele,frhs,afun,kappa,description)

           if varexist('frhs') this.problem_.frhs = frhs;
           else this.problem_.frhs = inline('[1]','x'); end;

           if varexist('afun') this.problem_.afun = afun;
           else this.problem_.afun = inline('[1*sign(x-0.5)+1]','x'); end;

           if varexist('kappa') this.problem_.kappa = kappa;
           else this.problem_.kappa = 0.00015625; end;

           if varexist('description') this.problem_.description = description;
           else this.problem_.description = 'Scatra 1d'; end;

           % gauss points
           this.gausspt_.nquad = 2;  % 2 gauss points
           this.gausspt_.wgt = [1; 1];
           this.gausspt_.xg = [-0.577350269; +0.577350269];

           % shape functions
           this.shpfct_.f1 = inline('-0.5*x+0.5','x');
           this.shpfct_.f2 = inline('0.5*x+0.5','x');
           this.shpfct_.f1d = inline('-0.5');
           this.shpfct_.f2d = inline('0.5');

           this.dis_.nElements = nele;

           numdofpernode = 1;
           numnodes = this.dis_.nElements + 1;
           numdofs = numnodes * numdofpernode;

           h = (this.xend_-this.xstart_)/this.dis_.nElements;

           xnod = zeros(numnodes,1);
           for i=1:numnodes
               xnod(i) = (i-1) * h;
           end

           this.problem_.dis = this.dis_;
           this.problem_.shpfct = this.shpfct_;
           this.problem_.gausspt = this.gausspt_;
           this.problem_.geo.h = h;
           this.problem_.geo.p = xnod;
           this.problem_.numnodes = numnodes;
           this.problem_.numdofpernode = numdofpernode;
           this.problem_.numdofs = numdofs;

           clear numdofpernode; clear numnodes; clear numdofs; clear h; clear xnod; clear i; clear nele;
       end

       function [A,b,problem] = Build(this)

           K = this.assemble_K();
           C = this.assemble_C();
           F = this.assemble_F();
           [Cstab,epe] = this.assemble_Cstab();
           Fstab = this.assemble_Fstab();

           A = K + C + Cstab;
           b = F + Fstab;

           A(1,:) = zeros(size(A,2),1);
           A(size(A,1),:) = zeros(size(A,2),1);
           A(1,1) = 1;     % set Dirichlet rb
           A(size(A,1),size(A,2)) = 1;
           b(1) = 0;
           b(length(b)) = 0;

           clear K; clear C; clear F; clear Fstab; clear Cstab;

           newMap = Map(size(A,1), 1); % generate full map (point-wise)
           A = Operator(A,newMap,newMap,@MatlabApply, this.problem_.description);

           clear newMap;

           problem = this.problem_;
       end

       function plot_vector(this,x)
           plot(this.problem_.geo.p,x);
       end

   end

   methods (Static = true)
   end

   methods (Access = private)
       function [K] = assemble_K(this)

           K = sparse(this.problem_.numdofs,this.problem_.numdofs);

           nel = this.problem_.dis.nElements;

           % loop over all elements
           for ele = 1:nel

               estif = zeros(2*this.problem_.numdofpernode);

               for iquad = 1:this.problem_.gausspt.nquad
                   % Gauss parameter * scalar jacobi determinant
                   fac =  this.problem_.gausspt.wgt(iquad) * this.problem_.geo.h / 2;

                   % global 1. derivatives of shape functions grad(phi)/grad(w) =
                   % dphi/dxi * 1/detJ
                   dn1dx = this.problem_.shpfct.f1d(0) * 2 / this.problem_.geo.h;
                   dn2dx = this.problem_.shpfct.f2d(0) * 2 / this.problem_.geo.h;

                   % shape functions at current gauss pt
                   n1 = this.problem_.shpfct.f1(this.problem_.gausspt.xg(iquad));
                   n2 = this.problem_.shpfct.f2(this.problem_.gausspt.xg(iquad));

                   %            estiv(1,1) = estif(1,1) + n1*problem.afun(x)*dn1dx*fac;
                   %            estiv(1,2) = estif(1,2) + n1*problem.afun(x)*dn2dx*fac;
                   %            estiv(2,1) = estif(2,1) + n2*problem.afun(x)*dn1dx*fac;
                   %            estiv(2,2) = estif(2,2) + n2*problem.afun(x)*dn2dx*fac;

                   estif(1,1) = estif(1,1) + dn1dx * this.problem_.kappa * dn1dx * fac;
                   estif(1,2) = estif(1,2) + dn1dx * this.problem_.kappa * dn2dx * fac;
                   estif(2,1) = estif(2,1) + dn2dx * this.problem_.kappa * dn1dx * fac;
                   estif(2,2) = estif(2,2) + dn2dx * this.problem_.kappa * dn2dx * fac;
               end

               % assembly
               lm = [ele,ele+1];
               K(lm(1),lm(1))=K(lm(1),lm(1))+estif(1,1);
               K(lm(1),lm(2))=K(lm(1),lm(2))+estif(1,2);
               K(lm(2),lm(1))=K(lm(2),lm(1))+estif(2,1);
               K(lm(2),lm(2))=K(lm(2),lm(2))+estif(2,2);
           end
       end % end assemble_K

       function [C] = assemble_C(this)

           C = sparse(this.problem_.numdofs,this.problem_.numdofs);

           nel = this.problem_.dis.nElements;

           % loop over all elements
           for ele = 1:nel

               estif = zeros(2*this.problem_.numdofpernode);

               for iquad = 1:this.problem_.gausspt.nquad
                   % Gauss parameter * scalar jacobi determinant
                   fac =  this.problem_.gausspt.wgt(iquad) * this.problem_.geo.h / 2;

                   % global 1. derivatives of shape functions grad(phi)/grad(w) =
                   % dphi/dxi * 1/detJ
                   dn1dx = this.problem_.shpfct.f1d(0) * 2 / this.problem_.geo.h;
                   dn2dx = this.problem_.shpfct.f2d(0) * 2 / this.problem_.geo.h;

                   % shape functions at current gauss pt
                   n1 = this.problem_.shpfct.f1(this.problem_.gausspt.xg(iquad));
                   n2 = this.problem_.shpfct.f2(this.problem_.gausspt.xg(iquad));

                   nodes = [ele,ele+1]; % nodes from 1d domain
                   coords = this.problem_.geo.p(nodes);

                   xgausp = (this.problem_.gausspt.xg(iquad)+1)*(this.problem_.geo.h/2) + this.problem_.geo.p(ele);
                   a = this.problem_.afun(xgausp);
                   estif(1,1) = estif(1,1) + n1*a*dn1dx*fac;
                   estif(1,2) = estif(1,2) + n1*a*dn2dx*fac;
                   estif(2,1) = estif(2,1) + n2*a*dn1dx*fac;
                   estif(2,2) = estif(2,2) + n2*a*dn2dx*fac;

               end

               % assembly
               lm = [ele,ele+1];
               C(lm(1),lm(1))=C(lm(1),lm(1))+estif(1,1);
               C(lm(1),lm(2))=C(lm(1),lm(2))+estif(1,2);
               C(lm(2),lm(1))=C(lm(2),lm(1))+estif(2,1);
               C(lm(2),lm(2))=C(lm(2),lm(2))+estif(2,2);
           end
       end % end assemble_C

       function [C,nele_epe] = assemble_Cstab(this)

           C = sparse(this.problem_.numdofs,this.problem_.numdofs);

           nel = this.problem_.dis.nElements;

           nele_epe = zeros(nel,1);    % vector with element peclet numbers

           % loop over all elements
           for ele = 1:nel

               estif = zeros(2*this.problem_.numdofpernode);

               for iquad = 1:this.problem_.gausspt.nquad
                   % Gauss parameter * scalar jacobi determinant
                   fac =  this.problem_.gausspt.wgt(iquad) * this.problem_.geo.h / 2;

                   % global 1. derivatives of shape functions grad(phi)/grad(w) =
                   % dphi/dxi * 1/detJ
                   dn1dx = this.problem_.shpfct.f1d(0) * 2 / this.problem_.geo.h;
                   dn2dx = this.problem_.shpfct.f2d(0) * 2 / this.problem_.geo.h;

                   % shape functions at current gauss pt
                   n1 = this.problem_.shpfct.f1(this.problem_.gausspt.xg(iquad));
                   n2 = this.problem_.shpfct.f2(this.problem_.gausspt.xg(iquad));

                   nodes = [ele,ele+1]; % nodes from 1d domain
                   coords = this.problem_.geo.p(nodes);

                   xgausp = (this.problem_.gausspt.xg(iquad)+1)*(this.problem_.geo.h/2) + this.problem_.geo.p(ele);

                   % calculate stabilization parameter
                   epe = this.problem_.afun(xgausp) * this.problem_.geo.h / (2*this.problem_.kappa);  % element peclet number
                   nele_epe(ele) = epe;
                   fun = coth(epe) - 1.0/epe;
                   %fun = min(1,epe/3);
                   if epe~= 0.0
                       tau = fun * this.problem_.geo.h/(2*this.problem_.afun(xgausp));
                   else
                       tau = 0;
                   end

                   % evaluate a at gauss point
                   a = this.problem_.afun(xgausp);

                   estif(1,1) = estif(1,1) + a*dn1dx*tau*a*dn1dx*fac;
                   estif(1,2) = estif(1,2) + a*dn1dx*tau*a*dn2dx*fac;
                   estif(2,1) = estif(2,1) + a*dn2dx*tau*a*dn1dx*fac;
                   estif(2,2) = estif(2,2) + a*dn2dx*tau*a*dn2dx*fac;

               end

               % assembly
               lm = [ele,ele+1];
               C(lm(1),lm(1))=C(lm(1),lm(1))+estif(1,1);
               C(lm(1),lm(2))=C(lm(1),lm(2))+estif(1,2);
               C(lm(2),lm(1))=C(lm(2),lm(1))+estif(2,1);
               C(lm(2),lm(2))=C(lm(2),lm(2))+estif(2,2);
           end
       end % end assemble_Cstab

       function [F] = assemble_F(this)

           F = sparse(this.problem_.numdofs,1);

           nel = this.problem_.dis.nElements;

           % loop over all elements
           for ele = 1:nel

               frhs = zeros(2*this.problem_.numdofpernode,1);

               for iquad = 1:this.problem_.gausspt.nquad
                   % Gauss parameter * scalar jacobi determinant
                   fac =  this.problem_.gausspt.wgt(iquad) * this.problem_.geo.h / 2;

                   % global 1. derivatives of shape functions grad(phi)/grad(w) =
                   % dphi/dxi * 1/detJ
                   dn1dx = this.problem_.shpfct.f1d(0) * 2 / this.problem_.geo.h;
                   dn2dx = this.problem_.shpfct.f2d(0) * 2 / this.problem_.geo.h;

                   % shape functions at current gauss pt
                   n1 = this.problem_.shpfct.f1(this.problem_.gausspt.xg(iquad));
                   n2 = this.problem_.shpfct.f2(this.problem_.gausspt.xg(iquad));

                   nodes = [ele,ele+1]; % nodes from 1d domain
                   coords = this.problem_.geo.p(nodes);

                   xgausp = (this.problem_.gausspt.xg(iquad)+1)*(this.problem_.geo.h/2) + this.problem_.geo.p(ele);

                   r = this.problem_.frhs(xgausp);

                   frhs(1) = frhs(1) + n1 * r * fac;
                   frhs(2) = frhs(2) + n2 * r * fac;
               end

               % assembly
               lm = [ele,ele+1];
               F(lm(1))=F(lm(1))+frhs(1);
               F(lm(2))=F(lm(2))+frhs(2);

           end
       end % end assemble_F

       function [F] = assemble_Fstab(this)

           F = sparse(this.problem_.numdofs,1);

           nel = this.problem_.dis.nElements;

           % loop over all elements
           for ele = 1:nel

               frhs = zeros(2*this.problem_.numdofpernode);

               for iquad = 1:this.problem_.gausspt.nquad
                   % Gauss parameter * scalar jacobi determinant
                   fac = this.problem_.gausspt.wgt(iquad) * this.problem_.geo.h / 2;

                   % global 1. derivatives of shape functions grad(phi)/grad(w) =
                   % dphi/dxi * 1/detJ
                   dn1dx = this.problem_.shpfct.f1d(0) * 2 / this.problem_.geo.h;
                   dn2dx = this.problem_.shpfct.f2d(0) * 2 / this.problem_.geo.h;

                   % shape functions at current gauss pt
                   n1 = this.problem_.shpfct.f1(this.problem_.gausspt.xg(iquad));
                   n2 = this.problem_.shpfct.f2(this.problem_.gausspt.xg(iquad));

                   nodes = [ele,ele+1]; % nodes from 1d domain
                   coords = this.problem_.geo.p(nodes);

                   xgausp = (this.problem_.gausspt.xg(iquad)+1)*(this.problem_.geo.h/2) + this.problem_.geo.p(ele);

                   % evaluate rhs at gauss point
                   r = this.problem_.frhs(xgausp);

                   % calculate stabilization parameter
                   epe = this.problem_.afun(xgausp) * this.problem_.geo.h / (2*this.problem_.kappa);  % element peclet number
                   fun = coth(epe) - 1.0/epe;
                   %fun = min(1,epe/3);

                   if epe ~= 0
                       tau = fun * this.problem_.geo.h/(2*this.problem_.afun(xgausp));
                   else
                       tau = 0;
                   end

                   % evaluate a at gauss point
                   a = this.problem_.afun(xgausp);

                   frhs(1) = frhs(1) + a * dn1dx * tau * r * fac;
                   frhs(2) = frhs(2) + a * dn2dx * tau * r * fac;
               end

               % assembly
               lm = [ele,ele+1];
               F(lm(1))=F(lm(1))+frhs(1);
               F(lm(2))=F(lm(2))+frhs(2);
           end
       end % assemble_Fstab
   end% end methods (private)

   methods (Access = protected)
       function Copy_(this,src,mc)
           [cmd,data,mc] = this.CopyCmd_(src,mc);
           eval(cmd);
       end
   end % protected methods
end