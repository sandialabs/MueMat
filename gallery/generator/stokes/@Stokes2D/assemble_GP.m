% note: nabla p \cdot q = -p \nabla \cdot q = - p \div q
function [GP] = assemble_GP(problem)
    GP = sparse(problem.numnodes * 2, problem.numnodes); % only velocity degrees!
    nel = problem.dis.nElements;

    %% loop over all elements
    for ele = 1:nel
       nodesperelement = 3; % fixed for triangle elements (3 nodes per element)

       estif = zeros(nodesperelement * 2,nodesperelement);  % element matrix

       %% loop over all gauss points
       for iquad = 1:problem.gausspt.nquad
            % get shape functions and derivatives for current gauss pt
            [funct,derxy,fac] = Stokes2D.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,problem.gausspt,problem.shpfct,problem.dis);

            % fill 6x3 elment pressure gradient matrix
            for ui = 1:nodesperelement
                fuipp  = ui;    % column
                for vi = 1:nodesperelement
                    fvi  = 2*vi - 1;
                    fvip = 2*vi;

                    estif(fvi,fuipp)=estif(fvi,fuipp) + ...
                        -fac * funct(ui) * derxy(vi,1);
                    estif(fvip,fuipp) = estif(fvip,fuipp) + ...
                        -fac * funct(ui) * derxy(vi,2);
                end
            end
       end

       %% assemble element GP matrix to GP
        % node numbers of element
       nodes = problem.geo.t(1:3,ele);

       for k=1:nodesperelement
           GP(2*nodes(k)-1,nodes(1)) = GP(2*nodes(k)-1,nodes(1)) + estif(2*k-1,1);
           GP(2*nodes(k)  ,nodes(1)) = GP(2*nodes(k)  ,nodes(1)) + estif(2*k  ,1);
           GP(2*nodes(k)-1,nodes(2)) = GP(2*nodes(k)-1,nodes(2)) + estif(2*k-1,2);
           GP(2*nodes(k)  ,nodes(2)) = GP(2*nodes(k)  ,nodes(2)) + estif(2*k  ,2);
           GP(2*nodes(k)-1,nodes(3)) = GP(2*nodes(k)-1,nodes(3)) + estif(2*k-1,3);
           GP(2*nodes(k)  ,nodes(3)) = GP(2*nodes(k)  ,nodes(3)) + estif(2*k  ,3);
       end
    end

    % dont forget dirichlet bcs!

end