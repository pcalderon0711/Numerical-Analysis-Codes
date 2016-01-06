function x = myBackSub_PGBC(A,b)

n=length(A);
x=zeros(n,1);
x(n)=b(n)/A(n,n);
for i = n-1:-1:1
    x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end

% Prepared by Pio Gabrielle B Calderon