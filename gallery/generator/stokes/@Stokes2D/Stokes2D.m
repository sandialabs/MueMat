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

classdef Stokes2D

    properties (Access = protected)
        problem;
    end

    methods (Access = public)
        function [this] = Stokes2D()
            % empty constructor
        end
    end

    methods (Access = public)

        function [A,rhs,sprob] = Build(this)
            sprob = this.problem;
            K  = this.assemble_K(this.problem);     % viscosity
            GP = this.assemble_GP(this.problem);    % pressure gradient
            DU = this.assemble_DU(this.problem);    % continuity equation
            Z  = this.assemble_Z(this.problem);     % pressure stabilization

            ST = [K,GP;-DU, -Z];   % adapt sign in continuity equation for symmetry

            [frhs,grhs] = this.assemble_F(this.problem); % RHS velocity part

            F = zeros(size(K,1)+size(DU,1),1);
            F(1:size(K,1)) = frhs;
            F(size(K,1)+1:size(F,1),1) = grhs;

            % Dirichlet BC velx
            bcnodes = this.problem.bc.dirichnodes_velx;

            for bc = 1:size(bcnodes,2)
                 ST((bcnodes(bc)-1)*2+1,:) = 0;
                 ST((bcnodes(bc)-1)*2+1,(bcnodes(bc)-1)*2+1) = 1;

                 F((bcnodes(bc))*2-1,1) = this.problem.bc.dirichvals_velx(bc);
            end

            % Dirichlet BC vely
            bcnodes = this.problem.bc.dirichnodes_vely;

            for bc = 1:size(bcnodes,2)
                 ST((bcnodes(bc)-1)*2+2,:) = 0;
                 ST((bcnodes(bc)-1)*2+2,(bcnodes(bc)-1)*2+2) = 1;

                 F((bcnodes(bc))*2,1)=this.problem.bc.dirichvals_vely(bc);
            end

            % Dirichlet BC for pressure?
            if isfield(this.problem.bc,'dirichnodes_pressure')
                bcnodes = this.problem.bc.dirichnodes_pressure;

                for bc = 1:size(bcnodes,2)
                    ST(size(K,1) + bcnodes(bc),:) = 0;
                    ST(size(K,1) + bcnodes(bc),size(K,1) + bcnodes(bc)) = 1;
                    F(size(K,1) + bcnodes(bc),1) = this.problem.bc.dirichvals_pressure(bc);
                end
            end

            [A,b] = this.transform_xyxypp_to_xypxyp(ST,F,this.problem);

            A = Operator(A,3,3);
            rhs = b;
        end

        function plot_solution(this,x)
            subplot(2,1,1);
            this.plot_stokes2d(this.problem.geo.p,this.problem.geo.e,this.problem.geo.t,'flowdata',[x(1:3:length(x)),x(2:3:length(x))]);
            subplot(2,1,2);
            this.plot_stokes2d(this.problem.geo.p,this.problem.geo.e,this.problem.geo.t,'zdata',x(3:3:length(x)));
        end

    end

    methods (Abstract = true)
        [problem] = defineproblem(this,p,e,t);
    end

    methods(Static = true)
        [p,r,t] = poi_mesh(g,n1,n2);
        [x,y]=pde_igeom(dl,bs,s);
        h=plot_stokes2d(p,e,t,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14);
        [uxy,tn,al2,al3]=tri2grid(p,t,u,tn,al2,al3);

        [DU] = assemble_DU(problem);
        [F,Fp] = assemble_F(problem);
        [GP] = assemble_GP(problem);
        [K] = assemble_K(problem);
        [Z] = assemble_Z(problem);
        [tau] = calc_tau(problem,cur_element);
        [funct,derxy,fac] = EvalShapeFuncAndDerivAtIntPoint(iquad,ele,gausspt,shpfct,dis);
    end

    methods(Access = private)


       function [A,b] = transform_xyxypp_to_xypxyp(this,ST,F,problem)
       % internal service function
            numdofpernode = problem.numdofpernode;
            veldofspernode = numdofpernode - 1;
            n = size(ST,1);
            nu = veldofspernode/numdofpernode * n;
            np = 1/numdofpernode * n;

            A2 = sparse(n,n);
            b = zeros(n,1);

            % reorder lines
            for i = 1:np
                A2((i-1)*numdofpernode+1,:) = ST((i-1)*veldofspernode+1,:);
                A2((i-1)*numdofpernode+2,:) = ST((i-1)*veldofspernode+2,:);
                A2((i-1)*numdofpernode+3,:) = ST(nu+i,:);
                b((i-1)*numdofpernode+1)   = F((i-1)*veldofspernode+1);
                b((i-1)*numdofpernode+2)   = F((i-1)*veldofspernode+2);
                b((i-1)*numdofpernode+3)   = F(nu+i);
            end

            % reorder columns
            for j = 1:np
               A(:,(j-1)*numdofpernode+1) = A2(:,(j-1)*veldofspernode+1);
               A(:,(j-1)*numdofpernode+2) = A2(:,(j-1)*veldofspernode+2);
               A(:,(j-1)*numdofpernode+3) = A2(:,nu+j);
            end

            clear A2; clear b2;
       end
    end
end