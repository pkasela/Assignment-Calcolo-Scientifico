function A = K1(n,h,a11)
%discretization with finite centered difference of the 
%second derivative for both x and y
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
% spdiags creates a sparse diag Matrix

