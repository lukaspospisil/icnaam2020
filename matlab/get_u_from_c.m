function [ u,x ] = get_u_from_c( c,xi_s,Nplot )
%GET_U_FROM_C Summary of this function goes here
%   Detailed explanation goes here

hh = (xi_s(2)-xi_s(1))/Nplot; 
x = xi_s(1):hh:xi_s(end);
u = zeros(size(x));

n = length(xi_s) - 1;

% goal: compute function values in u_plot
for i=1:n-1 % through elements
    idx_left = (i-1)*Nplot+1:i*Nplot;
    idx_right = i*Nplot+1:(i+1)*Nplot;
    
    xi_l = xi_s(i); % x_{i-1} ("left")
    xi = xi_s(i+1); % x_i ("center")
    xi_r = xi_s(i+2); % x_{i+1} ("right")
    
    u(idx_left) = u(idx_left) + c(2*i)*((x(idx_left) - xi).^3/(xi_l - xi)^2 - 2*(x(idx_left) - xi).^2/(xi_l - xi) + (x(idx_left) - xi));
    u(idx_right) = u(idx_right) + c(2*i)*((x(idx_right) - xi).^3/(xi_r - xi)^2 - 2*(x(idx_right) - xi).^2/(xi_r - xi) + (x(idx_right) - xi));

    u(idx_left) = u(idx_left) + c(2*i-1)*(1 + 2*((x(idx_left) - xi)/(xi_l - xi)).^3 - 3*((x(idx_left) - xi)/(xi_l - xi)).^2);
    u(idx_right) = u(idx_right) + c(2*i-1)*(1 + 2*((x(idx_right) - xi)/(xi_r - xi)).^3 - 3*((x(idx_right) - xi)/(xi_r - xi)).^2);
    
end


end

