function L=my_cholesky_PGBC(A)
% Input
%		A : positive definite matrix
% Output
%		L : lower triangular matrix s.t. L*L'=A
dim=length(A);
L=zeros(dim);
for k=1:dim
	L(k,k)=sqrt(A(k,k)-L(k,1:k-1)*L(k,1:k-1)');
	for i=k+1:dim
		L(i,k)=(A(i,k)-L(i,1:k-1)*L(k,1:k-1)')/L(k,k);
	end
end

%Prepared by PIO CALDERON