function y = myHorner(a, x0)

i = 1;
y = x0;
b= a(1);
while i < length(a)
    y=b*x0+a(i+1);
    b=y;
    i=i+1;
end

% prepared by PIO GABRIELLE B CALDERON