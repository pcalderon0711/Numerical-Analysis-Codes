function x3 = myMuller_PGBC(f, x0, x1, x2, maxiter)

i=3;
table=zeros(i,2);
table(1:3,1:2)=[x0 f(x0); x1 f(x1); x2 f(x2)];
for i = 4:maxiter+3
    h0=x0-x2;
    h1=x1-x2;
    h2=h0*h1*(h0-h1);
    fx0=f(x0);
    fx1=f(x1);
    fx2=f(x2);
    a=(h1*(fx0-fx2)-h0*(fx1-fx2))/h2;
    b=(((h0^2)*(fx1-fx2))-((h1^2)*(fx0-fx1)))/h2;
    D=sqrt(b^2-4*a*f(x2));
    R1=b+D;
    R2=b-D;
    if abs(b+D)>=abs(b-D)
        x3=x2-(2*fx2/R1);
    else
        x3=x2-(2*fx2/R2);
    end
    table(i,1:2)=[x3 f(x3)];
    x0=x1;
    x1=x2;
    x2=x3;
end
display(table);

%Prepared by Pio Gabrielle B. Calderon