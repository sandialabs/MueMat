function [F,Fp] = assemble_F(problem)

    nodesperelement = 3;    % nodes per element
    nel = problem.dis.nElements;

    % rechte Seite
    F = zeros(problem.numnodes * 2,1);  % velocity part
    Fp = zeros(problem.numnodes, 1); % pressure part (if not zero)

    for ele = 1:nel
       % nodes of current element
       nodes = problem.geo.t(1:3,ele);

       f = zeros(nodesperelement*2,1);  % only velocity part of RHS
       fp = zeros(nodesperelement,1);   % only pressure part of RHS

       % calculate stabilization factors tau
       tau = Stokes2D.calc_tau(problem,ele);

       % loop over all gauss points
       for iquad = 1:problem.gausspt.nquad

           [funct,derxy,fac] = Stokes2D.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,problem.gausspt,problem.shpfct,problem.dis);

           % coordinates of the nodes
           coords1 = problem.geo.p(:,nodes(1));
           coords2 = problem.geo.p(:,nodes(2));
           coords3 = problem.geo.p(:,nodes(3));
           xgausp = zeros(problem.numdofpernode,1);
           xgausp = coords1 + problem.gausspt.xg(iquad,1) * (coords2-coords1) + problem.gausspt.xg(iquad,2) * (coords3-coords1);

           r = zeros(nodesperelement*2,1);
           r = problem.frhs(xgausp(1),xgausp(2));

           % only velocity part of RHS
           for dofindex=1:2 % alle dofindizes
                for vi=1:nodesperelement  % alle Testfunktionen
                    f((vi-1)*2+dofindex,1) = f((vi-1)*2+dofindex,1) + fac*funct(vi) * r((vi-1)*problem.numdofpernode+dofindex,1);

                    % Konvektionsgeschwindigkeit für aktuellen Punkt
                    %coord = problem.geo.p(:,nodes(vi));
                    %a = zeros(1,2);
                    %a = problem.afun(coord(1),coord(2));

                    % Stabilisierungsparameter
                    % konvektiv taufac * conv(vi) * rhsint
                    %f((vi-1)*problem.numdofpernode+dofindex,1) = f((vi-1)*problem.numdofpernode+dofindex,1) + ...
                    %    tau(dofindex) * fac * (a(1)*derxy(vi,1) + a(2)*derxy(vi,2)) * r((vi-1)*problem.numdofpernode+dofindex,1);

                    % terms for pressure stabilization
                    f((vi-1)*2+dofindex,1) = f((vi-1)*2+dofindex,1) - tau * ( r((vi-1)*problem.numdofpernode+dofindex,1) * derxy(vi,dofindex));
                end
           end

           % only pressure part of RHS
            for vi=1:nodesperelement  % alle Testfunktionen
                fp(vi,1) = fp(vi,1) + fac*funct(vi) * r((vi-1)*problem.numdofpernode+2,1);
            end

       end

       for dofindex=1:2
           F((nodes(1)-1)*2+dofindex,1) = F((nodes(1)-1)*2+dofindex,1) + f(dofindex,1);
           F((nodes(2)-1)*2+dofindex,1) = F((nodes(2)-1)*2+dofindex,1) + f(1*2+dofindex,1);
           F((nodes(3)-1)*2+dofindex,1) = F((nodes(3)-1)*2+dofindex,1) + f(2*2+dofindex,1);
       end
       Fp(nodes(1),1) = Fp(nodes(1),1) + fp(1,1);
       Fp(nodes(2),1) = Fp(nodes(2),1) + fp(2,1);
       Fp(nodes(3),1) = Fp(nodes(3),1) + fp(3,1);
    end


end