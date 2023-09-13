% See also: VariableCoef1DFullyConstrained.m

%% Problem size
n = 81;

%% For n=81: 9 roots points at the middle of aggregates of size 9: 5 14 23 32 41 50 59 68 77

%% Rho coefficients
rho = ones(n+1,1);

P_noJump = VariableCoef1DFullyConstrained(n,rho);

%% Rho coefficients: add jump
epsi = 1e-5;

rho(27) = epsi;
rho(28) = epsi;
rho(29) = epsi;

P_Jump = VariableCoef1DFullyConstrained(n,rho);

axes1 = axes('Parent', figure, 'XTick', [5 14 23 32 41 50], 'DataAspectRatio', [1 4/50 1]); % xtick @ root points
%axes('DataAspectRatio', [1 3 1])
xlim(axes1,[3 52]);
ylim(axes1,[-0.1 1.1]);
box(axes1,'on');
hold(axes1,'all');

% With Jump: print basis func 2:5
plot1 = plot(P_Jump(:,2:5), 'LineWidth',2);
set(plot1(1),'Color',[0.749019622802734 0 0.749019622802734]);
set(plot1(2),'Color',[0 0 1]);
set(plot1(3),'Color',[1 0 0]);
%set(plot1(4),'Color',); use default

% Without: print basis func 3:4
plot2 = plot(P_noJump(:,3:4), '--', 'LineWidth',1);
set(plot2(1),'Color',[0 0 1]);
set(plot2(2),'Color',[1 0 0]);


