function [L,U]=my_crout_PGBC(A)
dim=length(A);
L=zeros(dim);
U=eye(dim);
L(1,1)=A(1,1);
for i=2:dim
    U(i-1,i)=A(i-1,i)/L(i-1,i-1);
    L(i,i-1)=A(i,i-1);
    L(i,i)=A(i,i)-U(i-1,i)*L(i,i-1);
end

% Prepared by PIO CALDERON