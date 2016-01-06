f = inline('sin(x)','x')
a = 0
b = 6.28
N = 100
method = 'T'
I = myCompositeIntegrator_PGBC(f,a,b,N,method)
method = 'S'
I = myCompositeIntegrator_PGBC(f,a,b,N,method)
