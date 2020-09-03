%
% The Dual Formulation of Discretized Beam Bending Problem with Sliding and Swivel Friction
% supplementary Matlab code for paper published on conference ICNAAM 2020
% Lukas Pospisil, Michaela Bobkova
% Department of Mathematics, Faculty of Civil Engineering, VSB-TU Ostrava, Czech Republic
% published under MIT Licence, 2020
%
% figure_solutions.m - the code generates figure in the paper which
% demonstrates the solution for different friction coefficients
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

% define sliding and swivel friction for several problems
problems{1}.g1 = 1e2;
problems{1}.g2 = 1e2;
problems{1}.legend = '$g_1=10^2\textrm{N},g_2=10^2\textrm{N}$';
problems{1}.color = [1.0,0,0];
problems{1}.markerstyle = '.';
problems{1}.linespec = '-';

problems{2}.g1 = 5e4;
problems{2}.g2 = 1e2;
problems{2}.legend = '$g_1=5\times10^4\textrm{N},g_2=10^2\textrm{N}$';
problems{2}.color = [0,0,1.0];
problems{2}.markerstyle = '.';
problems{2}.linespec = '--';

problems{3}.g1 = 1e2;
problems{3}.g2 = 5e4;
problems{3}.legend = '$g_1=10^2\textrm{N},g_2=5\times10^4\textrm{N}$';
problems{3}.color = [0,0.5,0.5];
problems{3}.markerstyle = '.';
problems{3}.linespec = ':';

problems{4}.g1 = 5e4;
problems{4}.g2 = 5e4;
problems{4}.legend = '$g_1=5\times10^4\textrm{N},g_2=5\times10^4\textrm{N}$';
problems{4}.color = [0.5,0.5,0];
problems{4}.markerstyle = '.';
problems{4}.linespec = '-.';


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

% something about symmetricity issues
A_hat = 0.5*(A_hat + A_hat');

for idx_problem=1:length(problems)
    g1 = problems{idx_problem}.g1;
    g2 = problems{idx_problem}.g2;
    
    g = [g1;g2];

    % use Matlab QP solver
    % quadprog is minimizing "0.5*X'*H*X + f'*X"
    % TODO: maybe play with some algorithm options?
    options = optimoptions('quadprog','Display','off','algorithm','interior-point-convex');
    problems{idx_problem}.lambda = quadprog(A_hat,-b_hat,[],[],[],[],-g,g,zeros(size(b_hat)),options);

    % recover primal solution
    problems{idx_problem}.c = A\(b - B'*problems{idx_problem}.lambda); % solve system instead of computing the inverse
end

%% plot solution
disp('- plot solution')

figure
hold on

y_min = -Inf*ones(1,length(problems));
y_max = Inf*ones(1,length(problems));

Nplot = 20; % density of plot (how many points to plot between xi_s)
mylegend_text = cell(size(problems));
mylegend_object = zeros(size(problems));

for idx_problem=1:length(problems)
    [u_plot,x_plot] = get_u_from_c(problems{idx_problem}.c,xi_s,Nplot);

    % y-axis limits
    y_min(idx_problem) = min(u_plot)-0.05;
    y_max(idx_problem) = max(u_plot)+0.05;

    mylegend_object(idx_problem) = plot(x_plot,u_plot,'b',...
                        'Color',problems{idx_problem}.color,...,...   
                        'LineStyle',problems{idx_problem}.linespec,...
                        'linewidth',2.0);
    mylegend_text{idx_problem} = problems{idx_problem}.legend;
    
    % plot FEM nodes
    for i=1:length(xi_s)
        [~,j] = min(abs(xi_s(i) - x_plot));
        plot(x_plot(j(1)), u_plot(j(1)),'b.',...
                        'Marker',problems{idx_problem}.markerstyle,...
                        'Color',problems{idx_problem}.color,'MarkerSize',20);
    end
    
end

% plot constrained point
plot([xi_s(x_hat_idx),xi_s(x_hat_idx)],[min(y_min),max(y_max)],'k--','linewidth',2.0)
    
hx = xlabel('$x$','Interpreter','latex','FontSize',16);
hy = ylabel('$u(x)$','Interpreter','latex','FontSize',16);

set(hx, 'FontSize', 16); 
set(hy, 'FontSize', 16); 

l = legend(mylegend_object,mylegend_text);
set(l, 'interpreter', 'latex', 'location', 'southwest', 'fontsize', 14);

axis([min(xi_s),max(xi_s),min(y_min),max(y_max)])
hold off



