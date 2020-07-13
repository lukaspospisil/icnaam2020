%
% The Dual Formulation of Discretized Beam Bending Problem with Sliding and Swivel Friction
% supplementary Matlab code for paper published on conference ICNAAM 2020
% Lukas Pospisil, Michaela Bobkova
% Department of Mathematics, Faculty of Civil Engineering, VSB-TU Ostrava, Czech Republic
% published under MIT Licence, 2020
%
% figure_basis.m - the code generates figure in the paper which
% demonstrates the basis functions
%

clear all

n = 6; % index of the last node (x_0,...,x_n) = number of elements
l = 2;

h = l/n; % size of intervals
xi_s = 0:h:l; % nodes
N = 2*(n-1); % dim Vh = number of unknowns


%% plot solution
disp('- plot solution')

Nplot = 20; % density of plot (how many points to plot between xi_s)
c_s = eye(N);

figure
for i=1:n-1
    for j=1:2
        subplot(2,n-1,(j-1)*(n-1)+i);
        hold on
        
        c = c_s(:,2*i-2+j);
        [u_plot,x_plot] = get_u_from_c(c,xi_s,Nplot);
        
        plot(x_plot,u_plot,'b','linewidth',2.0);
        
        % plot FEM nodes
        for k=1:length(xi_s)
            [~,jj] = min(abs(xi_s(k) - x_plot));
            plot(x_plot(jj(1)), u_plot(jj(1)),'b.','markersize',15);
        end
        
        hx = xlabel('$x$','Interpreter','latex');
        hy = ylabel(['$\varphi_{' num2str(2*i-2+j) '}(x)$'],'Interpreter','latex');
        
        set(hx, 'FontSize', 12);
        set(hy, 'FontSize', 12);
        
        if j==1
            axis([min(xi_s),max(xi_s),-0.15,1.05])
        end
        if j==2
            axis([min(xi_s),max(xi_s),-0.15,0.15])
        end
        hold off
    end
end


