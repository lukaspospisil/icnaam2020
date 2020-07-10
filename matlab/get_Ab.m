function [A,b,H,G] = get_Ab(n,h,E,J,f)
%GETA Summary of this function goes here
%   Detailed explanation goes here

H = [24/h^3, 0; 0, 8/h];
G = [-12/h^3, -6/h^2; 6/h^2, 2/h];
D3 = sparse(2:n-1,1:n-2,ones(1,n-2),n-1,n-1);
A = kron(speye(n-1),H) + kron(D3,G) + kron(D3',G');
A = E*J*A;

b = zeros(2*(n-1),1);
b(1:2:end) = f*h;


end

