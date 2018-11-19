% test

A = [1 2 3; 4 5 6; 7 8 9];
n = size(A,1);
m = 2;

vec = [1 2 3];

[K,H] = Arnoldi_orthogonalization(vec,A,m,n);