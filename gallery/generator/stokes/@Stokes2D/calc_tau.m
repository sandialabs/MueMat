% Stabilisierung nach Franca und Valentin 2000
function [tau] = calc_tau(problem,cur_element)

    % evaluate shape function
    e1 = 1/3;
    e2 = 1/3;
    funct = zeros(3,1);
    funct(1) = problem.shpfct.f1(e1,e2);
    funct(2) = problem.shpfct.f2(e1,e2);
    funct(3) = problem.shpfct.f3(e1,e2);

    % calculate are of triangle with formula of Heron
    xyze = problem.dis.e(cur_element).xyze;  % coords of current element (3x2)
    a = (xyze(1,1)-xyze(2,1))^2 + (xyze(1,2)-xyze(2,2))^2;
    b = (xyze(2,1)-xyze(3,1))^2 + (xyze(2,2)-xyze(3,2))^2;
    c = (xyze(3,1)-xyze(1,1))^2 + (xyze(3,2)-xyze(1,2))^2;
    area = 0.25 * sqrt(2*a*b+2*b*c+2*c*a-a^2-b^2-c^2);

    hk = sqrt(area);        % h-const for current element

    mk = 1/3;               % element type constant for tau

    % Geschwindigkeit (hier a) am Elementzentrum auswerten (also bei (1/3,1/3) der
    % Referenzzelle)
%     velint = zeros(2,1);    % 2 dim Geschwindigkeitsvektor am Elementzentrum
%
%     iel = 3;        % Knoten pro Element
%     for j=1:iel
%         nodes = problem.geo.t(1:3,cur_element);
%         coord = problem.geo.p(:,nodes(j));
%         a = zeros(1,2);
%         %a = problem.afun(coord(1),coord(2));
%         velint(1,1) = velint(1,1) + funct(j) * a(1);    % hier eigentlich evel(i+(2*j)) mit i=1
%         velint(2,1) = velint(2,1) + funct(j) * a(2);    % hier eigentlich evel(i+2(*j)) mit i=2
%     end
%
%     % Norm der Geschwindigkeit
%     vel_norm = norm(velint);

    % only stationary problems!

    % alle Skalaren abarbeiten
    %tau = zeros(problem.numdofpernode,1);   % fuer jeden Knoten tau's bestimmen (und zwar numdofpernode viele!)
    %for k=1:problem.numdofpernode
        %epe2 = mk * vel_norm * hk / problem.kappa(k);
        %xi2 = max(epe2,1.0);
    %    tau(k) = ((hk^2) * mk)/(4*problem.kappa(k)*1);
    %end
    tau = 0.12* ((hk^2) * mk)/(4*problem.kappa*1); % this is for the
    %artificial stokes problem
    %tau = ((hk^2) * mk)/(4*problem.kappa*1); % this is the original
end