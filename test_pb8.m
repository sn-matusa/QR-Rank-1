m=8;
n=6;
A=randn(m,n);
u=randn(m,1);
v=randn(n,1);
[Q1,R1]=qr(A);
[Q,R,A1]=QR_rang1(A,u,v) ;
norm(Q*R-A1)
[Q,R]=qrupdate(Q1,R1,u,v);    %o a doua testate folosind functia qrupdate
norm(Q*R-A1)