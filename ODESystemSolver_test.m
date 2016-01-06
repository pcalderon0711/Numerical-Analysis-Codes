f={@(t,y1,y2,y3) y2*y3, @(t,y1,y2,y3) -y1*y3, @(t,y1,y2,y3) -0.51*y1*y2};
ya=[0,1,1];
a=0;
b=12;
M=500;

mode="E";
[T,Y] = myODESystemSolver_PGBC(f,a,b,ya,M,mode);
plot(T,Y)

mode="R";
[T,Y] = myODESystemSolver_PGBC(f,a,b,ya,M,mode);
plot(T,Y)
