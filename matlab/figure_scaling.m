%
% The Dual Formulation of Discretized Beam Bending Problem with Sliding and Swivel Friction
% supplementary Matlab code for paper published on conference ICNAAM 2020
% Lukas Pospisil, Michaela Bobkova
% Department of Mathematics, Faculty of Civil Engineering, VSB-TU Ostrava, Czech Republic
% published under MIT Licence, 2020
%
% figure_scaling.m - the code generates figure in the paper, which
% demonstrates the scaling of the code for different problem dimensions
%

clear all

%% given parameters of the problem
E = 2.15e11; % Young's modulus [Nm-2]
l = 2; % length of beam [m]
x_hat = 1.3; % point obstacle with given friction [m]
f = -5e4; % load function [N]
v = 0.02; % cross section height [m]
s = 0.02; % cross section width [m]

% define sliding and swivel friction for several problems
g1 = 1e2;
g2 = 5e4;
g = [g1;g2];

% moment of inertia of the cross-section
% https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
J = s*v^3/12;

n_s = [10, 1e5:1e5:1e6]; % index of the last node (x_0,...,x_n) = number of elements

% time_all,
% time_primal_assembly,
% time_dual_assembly,
% time_dual_solve,
% time_primal_solve
times = zeros(5,length(n_s));
timers = zeros(1,length(n_s),'uint64');

N_rand = 100; % number of runs for every dimension
for n_idx = 1:length(n_s)
    n = n_s(n_idx);
    disp(['- solving n = ' num2str(n) ' (' num2str(n_idx) ' of ' num2str(length(n_s)) ')'])
    
    for n_rand=1:N_rand
        if mod(n_rand,10) == 0
            disp(['  - n_rand = ' num2str(n_rand)]);
        end
        
        timers(1) = tic;
        
        h = l/n; % size of intervals
        xi_s = 0:h:l; % nodes
        N = 2*(n-1); % dim Vh = number of unknowns
        
        % assemble objects of discretized problem
        timers(2) = tic;
        [A,b,H,G] = get_Ab(n,h,E,J,f);
        times(2,n_idx) = times(2,n_idx) + toc(timers(2));
        
        % find constrained index
        [~,j] = min(abs(xi_s-x_hat));
        x_hat_idx = j(1);
        
        % indexing issues (Matlab is indexing from 1)
        j = (x_hat_idx-1);
        
        % assemble objects in dual problem
        timers(3) = tic;
        B = sparse([1,2],[2*j-1,2*j],[1,1],2,N);
        
        BAinv = B/A; % solve system instead of computing the inverse
        A_hat = BAinv*B';
        b_hat = BAinv*b;
        
        % something about symmetricity issues
        A_hat = 0.5*(A_hat + A_hat');
        times(3,n_idx) = times(3,n_idx) + toc(timers(3));
        
        % use Matlab QP solver
        % quadprog is minimizing "0.5*X'*H*X + f'*X"
        % TODO: maybe play with some algorithm options?
        options = optimoptions('quadprog','Display','off','algorithm','interior-point-convex');
        timers(4) = tic;
        lambda = quadprog(A_hat,-b_hat,[],[],[],[],-g,g,zeros(size(b_hat)),options);
        times(4,n_idx) = times(4,n_idx) + toc(timers(4));
        
        % recover primal solution
        timers(5) = tic;
        c = A\(b - B'*lambda); % solve system instead of computing the inverse
        times(5,n_idx) = times(5,n_idx) + toc(timers(5));
        
        times(1,n_idx) = times(1,n_idx) + toc(timers(1));
        
    end
    
end

times = times(:,2:end)/N_rand; % compute average time
n_s = n_s(2:end); % ignore first dimension because of heat-up?

%% plot times

figure
hold on

mycolors{1} = [0,0,1];
mycolors{2} = [1,0,0];
mycolors{3} = [0,0.8,0];
mycolors{4} = [0.5,0.5,0];
mycolors{5} = [0.8,0,0.3];

mylinespecs = {'-','--','-',':','-.'};
mymarkers = {'o','x','s','^','>'};

for i=1:size(times,1)
    plot(n_s,times(i,:),'bo-',...
        'Color',mycolors{i},...
        'LineStyle',mylinespecs{i},...
        'Marker',mymarkers{i},...
        'MarkerSize',8,...
        'linewidth',2.0);
end

hx = xlabel('$n$','Interpreter','latex');
hy = ylabel('time [s]','Interpreter','latex');

set(hx, 'FontSize', 16);
set(hy, 'FontSize', 16);

l = legend('total time','primal assembly','dual assembly','dual solve', 'primal solve');
set(l, 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 14);

hold off



