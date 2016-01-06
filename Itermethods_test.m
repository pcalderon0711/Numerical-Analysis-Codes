format rat
%A = [2 -1 0; -1 2 -1; 0 -1 2]
%b = [1; 2; 3]
%w = 0
%x_0 = [0; 0; 0]
%tol = 1e-6
%n = 100
%mode = 'J'
%[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)
%mode = 'G'
%[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)
%mode = 'S'
%w = 1
%[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)

A = [4 2 0; 2 4 -1; 0 -1 4]
b = [5;3;6]
w = 0 
x_0 = [0;0;0]
n=100
tol=10**-4
mode = 'J'
[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)
mode = 'G'
[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)
mode = 'S'
[T,c,x_star, k, x_seq] = myitermethods_PGBC(A, b, w, x_0, tol, n, mode)