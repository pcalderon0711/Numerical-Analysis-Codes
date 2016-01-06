function [T, c, xstar, k, xk] = myitermethods_PGBC(A, b, w, x0, tol, maxiter, mode)

% Input
%		A  	    :  square coefficient matrix A s.t. Ax = b
%		b	    :  column matrix b s.t. Ax = b
%		w       :  relaxation factor for SOR method
%		x0      :  initial approximation for x written as a column vector
%		tol     :  error tolerance
%		maxiter :  maximum number of iterations to be performed
%		mode    :  'J' for Jacobi, 'G' for Gauss-Seidel, 'S' for SOR method
% Output
%		T       :  iteration matrix
%       c       :  column matrix for iteration
%		xstar   :  approximate solution to Ax = b
%		k       :  iteration number wherein x_star is obtained
%				:  the initial iterate x0 is not considered in the count
%		xk      :  sequence of iterates, first row is x0

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
switch mode
	case {'J'}
		T = -inv(D)*(L+U);
		c = inv(D)*b;
	case {'G'}
		T = -inv(D+L)*U;
		c = inv(D+L)*b;
	case {'S'}
		Tj = -inv(D)*(L+U);
		pj = max(abs(eig(Tj)));
		wOpt = 2/(1+sqrt(1-pj**2));
		if w != wOpt 
			fprintf('Input w is not optimal.');
			fprintf('Proceeding with the optimal w...\n');
			w = wOpt
		end
		T = inv(D+w*L)*((1-w)*D-w*U);
		c = inv(D+w*L)*w*b;
	otherwise
		fprintf('Please enter valid method: J or G or S');
		return
end

dim = length(x0);

% Initialize error and sequence of iterates
err = tol + 1;
% Transpose x0 since it is a column vector
% Row k+1 of xk corresponds to the kth iterate
xk = x0';

k = 1;
while err > tol && k < maxiter
	% xprev is previous iterate
	xprev = xk(k,1:dim)';
	x = T*xprev+c;
	% Error is calculated as the relative error wrt the previous iteration
	err = norm(x-xprev,inf)/norm(x);
	k = k + 1;
	xk(k,1:dim) = x';
end

xstar = xk(k,1:dim);
k = k-1;

% Prepared by Pio Calderon 2010-09736