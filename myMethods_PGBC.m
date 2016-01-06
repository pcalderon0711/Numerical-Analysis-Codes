function [k, xstar, abserr, relerr, xk] = myMethods_PGBC(mode, f,... 
        fprime, x0, x1, maxiter, tol)
    
% Input
%		mode	: method used to solve f(x) = 0, 'N' for 
%                 Newton's Method or 'S' for Secant Method
%		f       : inline function to be solved
%		fprime  : derivative of f as an inline function, ignored 
%                 for Secant Method
%		x0      : initial approximation for Newton's Method, first
%                 approximation for Secant Method
%		x1      : second approximation for Secant Method, ignored 
%                 for Newton's Method
%		maxiter : maximum number of iterations
%		tol     : tolerance for absolute error
% Outputs
% 		k		: the first iteration number wherein xk satisfies the
%                 stopping condition
%		xstar	: root approximation which satisfies the stopping 
%                 condition
%		abserr  : absolute error of xstar
%		relerr  : relative error of xstar
%		xk		: column matrix containing the sequence of root 
%                 approximations

i = 2;
% Initialize vector containing the sequence of approximations
xk = zeros(i,1);
% Initialize table which will be displayed at the end
table = zeros(i,4);

% Check if mode is Newton or Secant and check the validity of
% the inputs.
switch mode
	case {'N'}
		if f(x0) == 0
			fprintf('f(x0) = 0. x0 is a root of f.')
			return
		end
		if fprime(x0) == 0
			fprintf('fprime(x0) = 0. Please input new x0')
			return
		end
		xk(1) = x0;
        table(1,1:2) = [0,x0];
	case {'S'}
		if f(x0) == f(x1)
			fprintf('f(x0) = f(x1). Please input new x0 or x1.')
			return
		end
		if f(x0) == 0
			fprintf('f(x0) = 0. x0 is a root of f.')
			return
		elseif f(x1) == 0
			fprintf('f(x1) = 0. x1 is a root of f.')
			return
		end
		i = 3;
		xk(1:2) = [x0, x1];
        table(1:2, 1:2) = [0, x0; 1, x1];
		fx0 = f(x0);
		fx1 = f(x1);
    otherwise
		fprintf('Please enter valid mode: N or S.')
		return
end

while true
	switch mode
		case {'N'}
            % Assign new approx value to old approx variable.
			if i ~= 2
				x0 = x1;
            end
            % Compute new approximation.
			x1 = x0 - (f(x0)/fprime(x0));
			xk(i) = x1;
            % Compute errors
			abserr = abs(x1 - x0);
			relerr = abs((x1 - x0) / x1);
            % Fill up table with information.
            table(i, 1:4) = [i-1, x1, abserr, relerr];
		case {'S'}
            % Assign new approx values to old approx variables.
			if i ~= 3
				x0 = x1;
				x1 = x2;
				fx0 = fx1;
				fx1 = f(x2);
            end
            % Compute new approximation
			x2 = (x0*fx1 - x1*fx0) / (fx1 - fx0);
			xk(i) = x2;
            % Compute errors.
			abserr = abs(x2 - x1);
			relerr = abs((x2 - x1) / x2);
            % Fill up table with information.
            table(i, 1:4) = [i-1, x2, abserr, relerr];
    end
    % Check if current approximation satisfies the stopping conditions.
	if abserr < tol || relerr < tol || i == maxiter
        % Check if maxiter is reached and warn the user.
		if i == maxiter
			fprintf('Maximum number of iterations is reached.')
		end
		break;
	end
	i = i + 1;
end

% i-1 is the last iteration since list indexing starts with 1.
k = i - 1;
xstar = xk(i);

display(table)

% Prepared by Pio Gabrielle B. Calderon 2010-09736 Math 271.1 WF 4-5:30p