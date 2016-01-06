function [T,Y] = myODESystemSolver_PGBC(f,a,b,ya,M,mode)

% ODE System Solver
%
% Input
%			f    : a cell array of anonymous functions (@ notation)
%				 : containing the system of ODEs; see test function
%				 : for an example of usage
%			a    : the point at which the initial conditions are given
%			b    : the point up to which we solve the ODEs
%			ya   : the vector of initial conditions for the system f
%			M    : step size
%			mode : "E" for Euler and "R" for RK4
% Output
%			T    : vector of x-coordinates of solution
%			Y    : vector of y-coordinates of solution

numVar = length(ya);
h = (b-a)/M;
Y = zeros(M+1,numVar);
T = a:h:b;
Y(1,1:numVar) = ya;

if mode == 'E'
	for t=1:M
		for i=1:numVar
			% We convert Y(t,1:numVar) to a cell array with num2cell, then
			% we use the {:} operator to convert it to a comma-septd list
			Y(t+1,i)=Y(t,i)+h*f{i}(T(t),num2cell(Y(t,1:numVar)){:});
		end
	end
elseif mode == 'R'
	for t=1:M
		for i=1:numVar
			f1=f{i}(T(t),num2cell(Y(t,1:numVar)){:});
			f2=f{i}(T(t)+h/2,num2cell(Y(t,1:numVar)+(h*f1)/2){:});
			f3=f{i}(T(t)+h/2,num2cell(Y(t,1:numVar)+(h*f2)/2){:});
			f4=f{i}(T(t)+h,num2cell(Y(t,1:numVar)+h*f3){:});
			Y(t+1,i)=Y(t,i)+(h/6)*(f1+2*f2+2*f3+f4);
		end
	end
else
	fprintf("Please enter valid mode: E or R.")
end