function [K] = assemble_K(problem)
    K = sparse(problem.numnodes * 2, problem.numnodes * 2); % only velocity degrees!
    nel = problem.dis.nElements;

    %% loop over all elements
    for ele = 1:nel
       nodesperelement = 3; % fixed for triangle elements (3 nodes per element)
       kappa = problem.kappa;   % diffusivity constant
       estif = zeros(nodesperelement * 2);  % element matrix

       %% loop over all gauss points
       for iquad = 1:problem.gausspt.nquad
            % get shape functions and derivatives for current gauss pt
            [funct,derxy,fac] = Stokes2D.EvalShapeFuncAndDerivAtIntPoint(iquad,ele,problem.gausspt,problem.shpfct,problem.dis);

            % fill 6x6 elment viscosity matrix
            for ui = 1:nodesperelement
                fui  = 2*ui - 1;
                fuip = 2*ui;
                for vi = 1:nodesperelement
                    fvi  = 2*vi - 1;
                    fvip = 2*vi;

                    estif(fvi,fui)=estif(fvi,fui) + ...
                        fac * kappa * (2 * derxy(ui,1) * derxy(vi,1) + derxy(ui,2) * derxy(vi,2));
                    estif(fvi, fuip) = estif(fvi,fuip) + ...
                        fac * kappa * derxy(ui,1) * derxy(vi,2);
                    estif(fvip,fui) = estif(fvip,fui) + ...
                        fac * kappa * derxy(vi,1) * derxy(ui,2);
                    estif(fvip,fuip) = estif(fvip,fuip) + ...
                        fac * kappa * (derxy(ui,1)*derxy(vi,1) + 2*derxy(ui,2)*derxy(vi,2));
                end
            end
       end

       %% assemble element viscosity matrix in K
       % node numbers of element
       nodes = problem.geo.t(1:3,ele);
       % node 1
       K(2*nodes(1)-1,2*nodes(1)-1) = K(2*nodes(1)-1,2*nodes(1)-1) + estif(1,1);
       K(2*nodes(1)-1,2*nodes(1)  ) = K(2*nodes(1)-1,2*nodes(1)  ) + estif(1,2);
       K(2*nodes(1)  ,2*nodes(1)-1) = K(2*nodes(1)  ,2*nodes(1)-1) + estif(2,1);
       K(2*nodes(1)  ,2*nodes(1)  ) = K(2*nodes(1)  ,2*nodes(1)  ) + estif(2,2);

       K(2*nodes(1)-1,2*nodes(2)-1) = K(2*nodes(1)-1,2*nodes(2)-1) + estif(1,3);
       K(2*nodes(1)-1,2*nodes(2)  ) = K(2*nodes(1)-1,2*nodes(2)  ) + estif(1,4);
       K(2*nodes(1)  ,2*nodes(2)-1) = K(2*nodes(1)  ,2*nodes(2)-1) + estif(2,3);
       K(2*nodes(1)  ,2*nodes(2)  ) = K(2*nodes(1)  ,2*nodes(2)  ) + estif(2,4);

       K(2*nodes(1)-1,2*nodes(3)-1) = K(2*nodes(1)-1,2*nodes(3)-1) + estif(1,5);
       K(2*nodes(1)-1,2*nodes(3)  ) = K(2*nodes(1)-1,2*nodes(3)  ) + estif(1,6);
       K(2*nodes(1)  ,2*nodes(3)-1) = K(2*nodes(1)  ,2*nodes(3)-1) + estif(2,5);
       K(2*nodes(1)  ,2*nodes(3)  ) = K(2*nodes(1)  ,2*nodes(3)  ) + estif(2,6);

       % node 2
       K(2*nodes(2)-1,2*nodes(1)-1) = K(2*nodes(2)-1,2*nodes(1)-1) + estif(3,1);
       K(2*nodes(2)-1,2*nodes(1)  ) = K(2*nodes(2)-1,2*nodes(1)  ) + estif(3,2);
       K(2*nodes(2)  ,2*nodes(1)-1) = K(2*nodes(2)  ,2*nodes(1)-1) + estif(4,1);
       K(2*nodes(2)  ,2*nodes(1)  ) = K(2*nodes(2)  ,2*nodes(1)  ) + estif(4,2);

       K(2*nodes(2)-1,2*nodes(2)-1) = K(2*nodes(2)-1,2*nodes(2)-1) + estif(3,3);
       K(2*nodes(2)-1,2*nodes(2)  ) = K(2*nodes(2)-1,2*nodes(2)  ) + estif(3,4);
       K(2*nodes(2)  ,2*nodes(2)-1) = K(2*nodes(2)  ,2*nodes(2)-1) + estif(4,3);
       K(2*nodes(2)  ,2*nodes(2)  ) = K(2*nodes(2)  ,2*nodes(2)  ) + estif(4,4);

       K(2*nodes(2)-1,2*nodes(3)-1) = K(2*nodes(2)-1,2*nodes(3)-1) + estif(3,5);
       K(2*nodes(2)-1,2*nodes(3)  ) = K(2*nodes(2)-1,2*nodes(3)  ) + estif(3,6);
       K(2*nodes(2)  ,2*nodes(3)-1) = K(2*nodes(2)  ,2*nodes(3)-1) + estif(4,5);
       K(2*nodes(2)  ,2*nodes(3)  ) = K(2*nodes(2)  ,2*nodes(3)  ) + estif(4,6);

       % node 3
       K(2*nodes(3)-1,2*nodes(1)-1) = K(2*nodes(3)-1,2*nodes(1)-1) + estif(5,1);
       K(2*nodes(3)-1,2*nodes(1)  ) = K(2*nodes(3)-1,2*nodes(1)  ) + estif(5,2);
       K(2*nodes(3)  ,2*nodes(1)-1) = K(2*nodes(3)  ,2*nodes(1)-1) + estif(6,1);
       K(2*nodes(3)  ,2*nodes(1)  ) = K(2*nodes(3)  ,2*nodes(1)  ) + estif(6,2);

       K(2*nodes(3)-1,2*nodes(2)-1) = K(2*nodes(3)-1,2*nodes(2)-1) + estif(5,3);
       K(2*nodes(3)-1,2*nodes(2)  ) = K(2*nodes(3)-1,2*nodes(2)  ) + estif(5,4);
       K(2*nodes(3)  ,2*nodes(2)-1) = K(2*nodes(3)  ,2*nodes(2)-1) + estif(6,3);
       K(2*nodes(3)  ,2*nodes(2)  ) = K(2*nodes(3)  ,2*nodes(2)  ) + estif(6,4);

       K(2*nodes(3)-1,2*nodes(3)-1) = K(2*nodes(3)-1,2*nodes(3)-1) + estif(5,5);
       K(2*nodes(3)-1,2*nodes(3)  ) = K(2*nodes(3)-1,2*nodes(3)  ) + estif(5,6);
       K(2*nodes(3)  ,2*nodes(3)-1) = K(2*nodes(3)  ,2*nodes(3)-1) + estif(6,5);
       K(2*nodes(3)  ,2*nodes(3)  ) = K(2*nodes(3)  ,2*nodes(3)  ) + estif(6,6);


    end

    % dont forget dirichlet bcs!

    % Dirichlet RB
%     bcnodes = problem.bc.dirichnodes;
%
%     for bc = 1 : size(bcnodes,2)
%         for dofindex=1:2 % alle dofindizes
%             K((bcnodes(bc)-1)*2+dofindex,:) = 0;
%             K((bcnodes(bc)-1)*2+dofindex,(bcnodes(bc)-1)*2+dofindex) = 1;
%         end
%     end

end