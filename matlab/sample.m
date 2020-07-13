%
% The Dual Formulation of Discretized Beam Bending Problem with Sliding and Swivel Friction
% supplementary Matlab code for paper published on conference ICNAAM 2020
% Lukas Pospisil, Michaela Bobkova
% Department of Mathematics, Faculty of Civil Engineering, VSB-TU Ostrava, Czech Republic
% published under MIT Licence, 2020
%
% sample.m - the code demonstrates how to solve one problem with given parameters
%

clear all

%% given parameters of the problem
E = 2.15e11; % Young's modulus [Nm-2]
l = 2; % length of beam [m]
x_hat = 1.3; % point obstacle with given friction [m]
f = -5e4; % load function [N]
v = 0.02; % cross section height [m]
s = 0.02; % cross section width [m]
n = 20; % index of the last node (x_0,...,x_n) = number of elements

g1 = 5e4; % sliding friction [N]
g2 = 1e2; % swivel friction [N]

% moment of inertia of the cross-section
% https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
J = s*v^3/12;

h = l/n; % size of intervals
xi_s = 0:h:l; % nodes
N = 2*(n-1); % dim Vh = number of unknowns

%% assemble objects of discretized problem
disp('- assembly objects in QP')
[A,b,H,G] = get_Ab(n,h,E,J,f);

% find constrained index
[~,j] = min(abs(xi_s-x_hat));
x_hat_idx = j(1);

%% compute coeffs in the basis
disp('- QP solution')

% indexing issues (Matlab is indexing from 1)
j = (x_hat_idx-1);

% assemble objects in dual problem
B = sparse([1,2],[2*j-1,2*j],[1,1],2,N);

BAinv = B/A; % solve system instead of computing the inverse
A_hat = BAinv*B';
b_hat = BAinv*b;
g = [g1;g2];

% something about symmetricity issues
A_hat = 0.5*(A_hat + A_hat');

% use Matlab QP solver
% quadprog is minimizing "0.5*X'*H*X + f'*X"
% TODO: maybe play with some algorithm options?
options = optimoptions('quadprog','Display','iter','algorithm','interior-point-convex');
lambda = quadprog(A_hat,-b_hat,[],[],[],[],-g,g,zeros(size(b_hat)),options);

% recover primal solution
c = A\(b - B'*lambda); % solve system instead of computing the inverse

%% plot solution
disp('- plot solution')

Nplot = 20; % density of plot (how many points to plot between xi_s)
[u_plot,x_plot] = get_u_from_c(c,xi_s,Nplot);

% y-axis limits
y_min = min(u_plot)-0.05;
y_max = max(u_plot)+0.05;

figure
hold on
plot(x_plot,u_plot,'b','linewidth',2.0)

% plot constrained point
plot([xi_s(x_hat_idx),xi_s(x_hat_idx)],[y_min,y_max],'k--','linewidth',2.0)

% plot FEM nodes
for i=1:length(xi_s)
    [~,j] = min(abs(xi_s(i) - x_plot));
    plot(x_plot(j(1)), u_plot(j(1)),'b.','markersize',15);
end

xlabel('$x$','Interpreter','latex','fontsize',12)
ylabel('$u(x)$','Interpreter','latex','fontsize',12)

axis([min(xi_s),max(xi_s),y_min,y_max])
hold off
