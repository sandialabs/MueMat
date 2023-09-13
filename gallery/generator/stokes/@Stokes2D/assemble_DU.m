function [DU] = assemble_DU(problem)
    DU = sparse(problem.numnodes, 2 * problem.numnodes);
    nel = problem.dis.nElements;

    %% loop over all elements
    for ele = 1:nel
       nodesperelement = 3; % fixed for triangle elements (3 nodes per element)

       estif = zeros(nodesperelement, 2*nodesperelement);  % element matrix

       %% loop over all gauss points
       for iquad = 1:problem.gausspt.nquad
            % get shape functions and derivatives for current gauss pt
            [funct,derxy,fac] = Stokes2D.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,problem.gausspt,problem.shpfct,problem.dis);

            % fill 3x6 elment divergence of U matrix
            for vi = 1:nodesperelement
                fvipp  = vi;    % column
                for ui = 1:nodesperelement
                    fui  = 2*ui - 1;
                    fuip = 2*ui;

                    estif(fvipp,fui)=estif(fvipp,fui) + ...
                        +fac * funct(vi) * derxy(ui,1);
                    estif(fvipp,fuip) = estif(fvipp,fuip) + ...
                        +fac * funct(vi) * derxy(ui,2);
                end
            end
       end

       %% assemble element GP matrix to GP
        % node numbers of element
       nodes = problem.geo.t(1:3,ele);

       for k=1:nodesperelement
           DU(nodes(1),2*nodes(k)-1) = DU(nodes(1),2*nodes(k)-1) + estif(1,2*k-1);
           DU(nodes(1),2*nodes(k)  ) = DU(nodes(1),2*nodes(k)  ) + estif(1,2*k  );
           DU(nodes(2),2*nodes(k)-1) = DU(nodes(2),2*nodes(k)-1) + estif(2,2*k-1);
           DU(nodes(2),2*nodes(k)  ) = DU(nodes(2),2*nodes(k)  ) + estif(2,2*k  );
           DU(nodes(3),2*nodes(k)-1) = DU(nodes(3),2*nodes(k)-1) + estif(3,2*k-1);
           DU(nodes(3),2*nodes(k)  ) = DU(nodes(3),2*nodes(k)  ) + estif(3,2*k  );
       end

    end

    % dont forget dirichlet bcs!


end