function [z,C,D] = myCubicSplines_PGBC(x,y)

%	Natural Cubic Splines Coefficient Calculator
%	
%	This program calculates the coefficients z_i, C_i, D_i 
%	described in the handout for the natural cubic splines
%	method.

%	Input
%		x : (1 x n) column vector of interpolation nodes x_i
%		y : (1 x n) column vector of f(x_i) in the same order as x
%	Output
%		z : (1 x n) column vector containing z_i 
%		C : (1 x n-1) column vector containing C_i
%		D : (1 x n-1) column vector containing D_i

numNodes = length(x);
z = zeros(numNodes,1);
h = zeros(numNodes-1,1);
h(1) = x(2)-x(1);

% Obtain A and b to form the system Az = b, which considers the
% continuity of the derivative of the splines at the internal nodes.
A = zeros(numNodes-2,numNodes-2);
dimA = length(A);
b = zeros(numNodes-2,1);
for i = 2:numNodes-1	
	h(i) = x(i+1) - x(i);
	% Here, we fill up the interior rows.
	if i > 2 && i < numNodes-1
		A(i-1,i-2:i) = [h(i-1), 2*(h(i-1) + h(i)), h(i)];
	% This considers the case when A has size 1.
	elseif i == 2 && i == numNodes-1
		A(i-1,i-1) = 2*(h(i-1) + h(i));
	% We do not include the coeff of z_0 in the matrix.
	elseif i == 2
		A(i-1,i-1:i) = [2*(h(i-1) + h(i)), h(i)];
	% We do not include the coeff of z_n in the matrix.
	elseif i == numNodes-1
		A(i-1,i-2:i-1) = [h(i-1), 2*(h(i-1) + h(i))];
	end
	b(i-1) = 6*(((y(i+1)-y(i))/h(i))-((y(i)-y(i-1))/h(i-1)));

end

% Solve for z_i, C_i and D_i using the formula in the handout.
z(2:numNodes-1) = A\b;
C = (1./h) .* (y(2:numNodes) - ((h.*h)/6).*z(2:numNodes));
D = (1./h) .* (y(1:numNodes-1) - ((h.*h)/6).*z(1:numNodes-1));

% Prepared by PIO CALDERON 2010-09736