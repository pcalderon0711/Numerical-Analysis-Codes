function I = myCompositeIntegrator_PGBC(f,a,b,N,method)

% Composite Trapezoidal/ Simpson Integrator
% Integrate the given function over the interval [a,b].

%
% Input
%		f      : inline function to be integrated
%		a      : left endpoint of interval [a,b]
%		b      : right endpoint of interval [a,b]
%		N      : number of sub-intervals to divide [a,b] into
%		method : 'T' for composite Trapezoidal and 'S' for composite Simpson
% Output
%		I      : approximate value of the integral of f over [a,b]

h = (b-a)/N;
x = a:h:b;
I = 0;

if method == 'T'
	for i = 1:N-1
		I = I + h*f(x(i));
	end
	I = I + (h/2) * (f(a) + f(b));
elseif method == 'S'
	for i = 1:N-1
		I = I + f(x(i)) + 4*f((x(i)+x(i+1))/2) + f(x(i+1));
	end
	I = I*(h/6);
else
	fprintf('Enter valid method: T or S.')
	return;
end

% PREPARED BY PIO CALDERON