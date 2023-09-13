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

classdef Stokes2D_drivencavity < Stokes2D

    methods(Access = public)

        function [this] = Stokes2D_drivencavity()
            load geometry.mat
            [p,e,t] = Stokes2D.poi_mesh (g,10,10);

            this.problem = this.defineproblem(p,e,t);

        end

        function [problem] = defineproblem(this,p,e,t)

            % define gauss points
            gausspt.nquad = 3;
            gausspt.wgt = [1/6; 1/6; 1/6];
            gausspt.xg = [1/6 1/6; 2/3 1/6; 1/6 2/3];

            % define shape functions
            shpfct.f1 = inline('1-r-s','r','s');
            shpfct.f2 = inline('r','r','s');
            shpfct.f3 = inline('s','r','s');
            shpfct.f1d= inline('[-1;-1]','r','s');
            shpfct.f2d= inline('[1;0]','r','s');
            shpfct.f3d= inline('[0;1]','r','s');

            % number of elements (number of triangles)
            dis.nElements = size(t,2);
            for iele=1:dis.nElements
               temp = p(:,t(1:3,iele));
               ele.xyze = zeros(3,2);
               ele.xyze(1,:) = temp(:,1)';  % coordinates of node 1 in element
               ele.xyze(2,:) = temp(:,2)';  % coordinates of node 2 in element
               ele.xyze(3,:) = temp(:,3)';  % coordinates of node 3 in element
               dis.e(iele) = ele;
               clear temp;
               clear ele;
            end
            clear iele;

            % number of dofs per node
            % 2 velocity dofs and 1 pressure dof
            numdofpernode = 3;

            % number of nodes and dofs
            numnodes = size(p(1,:),2);
            numdofs = numnodes * numdofpernode;

            problem.dis = dis;
            problem.shpfct = shpfct;
            problem.gausspt = gausspt;
            problem.geo.p = p;
            problem.geo.e = e;
            problem.geo.t = t;
            problem.numnodes = numnodes;
            problem.numdofpernode = numdofpernode;
            problem.numdofs = numdofs;

            % 1d example
            problem.frhs = inline('[1;0;1;1;0;1;1;0;1]','x','y');  % 3 nodes (a la numdofpernode
            %dofs!) per element!

            % zero boundaries
            problem.bc.dirichnodes_velx = find(p(1,:)==0 | p(2,:)==0 | p(1,:)==1);% | p(2,:)==1);
            problem.bc.dirichvals_velx = zeros(length(problem.bc.dirichnodes_velx),1);
            problem.bc.dirichnodes_vely = find(p(1,:)==0 | p(2,:)==0 | p(1,:)==1);% | p(2,:)==1);
            problem.bc.dirichvals_vely = zeros(length(problem.bc.dirichnodes_vely),1);

            % boundary of driven lid cavity
            temp = find(p(2,:)==1);
            problem.bc.dirichnodes_velx = [problem.bc.dirichnodes_velx,temp];
            problem.bc.dirichnodes_vely = [problem.bc.dirichnodes_vely,temp];
            tempval = zeros(2*length(temp),1);

            for k=1:length(temp)
                % current top boundary node
                nodeidx = temp(k);
                xc = p(1,nodeidx);
                %tempval(1:2:length(tempval))=1;
                tempval(2*k-1) = xc^2 - 2*xc^3 + xc^4;
            end
            clear nodeidx; clear xc;

            problem.bc.dirichvals_velx = [problem.bc.dirichvals_velx;tempval];
            problem.bc.dirichvals_vely = [problem.bc.dirichvals_vely;zeros(length(tempval),1)];

            problem.bc.dirichnodes_pressure = 13;
            problem.bc.dirichvals_pressure = 0.0;

            clear temp; clear tempval;

            % diffusivity
            problem.kappa = 1; %ones(problem.numdofpernode,1);

            problem.description = 'Stokes problem';
        end
    end

end