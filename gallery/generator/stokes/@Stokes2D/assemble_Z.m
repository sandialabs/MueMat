function [Z] = assemble_Z(problem)
    Z = sparse(problem.numnodes, problem.numnodes); % only pressure dofs
    nel = problem.dis.nElements;

    %% loop over all elements
    for ele = 1:nel
       nodesperelement = 3; % fixed for triangle elements (3 nodes per element)

       estif = zeros(nodesperelement);  % element matrix

       % calculate tau for stabilization of pressure
       tau = Stokes2D.calc_tau(problem,ele);

       %% loop over all gauss points
       for iquad = 1:problem.gausspt.nquad
            % get shape functions and derivatives for current gauss pt
            [funct,derxy,fac] = Stokes2D.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,problem.gausspt,problem.shpfct,problem.dis);

            % fill 6x6 elment viscosity matrix
            for ui = 1:nodesperelement
                fuipp = ui;
                for vi = 1:nodesperelement
                    fvipp = vi;

                    estif(fvipp,fuipp)=estif(fvipp,fuipp) + ...
                        fac * tau * (derxy(ui,1) * derxy(vi,1) + derxy(ui,2) * derxy(vi,2));
                end
            end
       end

       %% assemble
       % node numbers of element
       nodes = problem.geo.t(1:3,ele);

       for i=1:nodesperelement
           for j=1:nodesperelement
               Z(nodes(i),nodes(j)) = Z(nodes(i),nodes(j)) + estif(i,j);
           end
       end

    end

end