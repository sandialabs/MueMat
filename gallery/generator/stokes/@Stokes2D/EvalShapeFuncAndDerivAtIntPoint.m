% current gauss point: iquad (1-3)
% current element: ele
% derxy: global derivative (3x2 matrix)
function [funct,derxy,fac] = EvalShapeFuncAndDerivAtIntPoint(iquad,ele,gausspt,shpfct,dis)
    % coordinates of current gauss pt
    e1 = gausspt.xg(iquad,1);
    e2 = gausspt.xg(iquad,2);

    % eval the three shape functions at current gauss pot (e1,e2)
    funct = zeros(3,1);
    funct(1) = shpfct.f1(e1,e2);
    funct(2) = shpfct.f2(e1,e2);
    funct(3) = shpfct.f3(e1,e2);

    % calculate derivatives of the shape functions at gauss pt (e1, e2)
    % -> (3x2 matrix)
    deriv = zeros(3,2);
    deriv(1,:) = shpfct.f1d(e1,e2)';
    deriv(2,:) = shpfct.f2d(e1,e2)';
    deriv(3,:) = shpfct.f3d(e1,e2)';

    iel = 3;    % number of nodes per element

    xjm = zeros(2,2);
    %      get Jacobian matrix and determinant
    %      actually compute its transpose....
    %            +-       -+ T      +-       -+
    %            | dx   dx |        | dx   dy |
    %            | --   -- |        | --   -- |
    %            | dr   ds |        | dr   dr |
    %     x_jm = |         |   =    |         |
    %            | dy   dy |        | dx   dy |
    %            | --   -- |        | --   -- |
    %            | dr   ds |        | ds   ds |
    %            +-       -+        +-       -+
    for i=1:2
        for j=1:2
            dum = 0;
            for l=1:iel
                dum = dum + deriv(l,i)*dis.e(ele).xyze(l,j);
            end
            xjm(i,j) = dum;
        end
    end

    % determinant of xjm with Sarrus (~> det (J) )
    determinante = xjm(1,1)*xjm(2,2) - xjm(1,2)*xjm(2,1);

    % check det < 0?
    if determinante < 0
        error('det < 0');
    end

    % weights of gauss integration * det(J) = factor for gauss integration
    fac = gausspt.wgt(iquad)*determinante;

    % inverse of xjm
    xij = zeros(2,2);
    xij(1,1) = +xjm(2,2)/determinante;
    xij(2,1) = -xjm(2,1)/determinante;
    xij(1,2) = -xjm(1,2)/determinante;
    xij(2,2) = +xjm(1,1)/determinante;

    %   --------------------------------------------------------------
    %        compute global first derivates using inverse xij
    %   --------------------------------------------------------------
    %    Use the Jacobian and the known derivatives in element coordinate
    %    directions on the right hand side to compute the derivatives in
    %    global coordinate directions
    %       +-          -+     +-    -+      +-    -+
    %       |  dx    dy  |     | dN_k |      | dN_k |
    %       |  --    --  |     | ---- |      | ---- |
    %       |  dr    dr  |     |  dx  |      |  dr  |
    %       |            |  *  |      |   =  |      | for all k
    %       |  dx    dy  |     | dN_k |      | dN_k |
    %       |  --    --  |     | ---- |      | ---- |
    %       |  ds    ds  |     |  dy  |      |  ds  |
    %       +-          -+     +-    -+      +-    -+
    %     `-------v-------´
    %            xjm             derxy         deriv
    %
    % global derivatives
    derxy = zeros(3,2);
    for k=1:iel
        derxy(k,1) = xij(1,1)*deriv(k,1) + xij(1,2)*deriv(k,2);
        derxy(k,2) = xij(2,1)*deriv(k,1) + xij(2,2)*deriv(k,2);
    end

end