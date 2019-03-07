% Ali Heydari
% Iterative methods
% GMRES Method for 9 point Laplacian w/ Dirchlet BC

f = inline('-5 * cos(x + 2 * y)','x','y');
g = inline('cos(x + 2*y)','x','y');


h = 0.1;
n = 1/h;
N = (n-1)^2;
m = n;
x_0 = zeros(N,1);

% initializing


A = zeros(N,N);
b = zeros(N,1);



[XM,YM] = meshgrid(x,y);
G = [];

for i = 1 : n

for j = 1 : n
G(i,j) = g(XM(i,j),YM(i,j));

end

end

Q = [];

for i = 1 : n

for j = 1 : n
Q(i,j) = f(XM(i,j),YM(i,j));

end

end



for j=1:n-1

yj = j*h;

%    make matrix A for the exact solution later on

    for i=1:n-1
        xi = i*h;

        k = i + (j-1)*(n-1);

        A(k,k) = -4;
            if i > 1
                A(k,k-1) = 1;
            end
            if i < n-1
                A(k,k+1) = 1 ;
            end
            if j > 1
                A(k,k-(n-1)) = 1;
            end

            if j < n-1
                A(k,k+(n-1)) = 1;
            end


    
    b(k) = h^2 * ( Q(i+1,j) + Q(i-1,j) + Q(i,j+1) + ...
                Q(i,j-1) + 8 * Q(i,j));
 
        % on the left
        if i==1
            b(k) = h^2 * b(k) - g(0,yj);
        end
        % on the right
        if i==n-1
            b(k) = h^2 * b(k) - g(1,yj);
        end
        % on the top
        if j==n-1
            b(k) = h^2 * b(k) - g(xi,1);
        end
        % on the bottom
        if j==1
            b(k) = h^2 * b(k) - g(xi,0);
        end
  
  end
  
  
end


    tol = 10^-7;
 [X,K,H,r_new,time,counter] = GMRES_1(x_0,A,b,m,tol);
   
    counter = counter * 2;
fprintf("iter: %i    time:%i\n" ,counter,time);




u = X;

%  Plot results.
figure(1);
ugrid = reshape(u,n-1,n-1);  % This makes the Nx1 vector u into an n-1 x n-1
                             % array.
surf([h:h:1-h]',[h:h:1-h]',ugrid)  %This plots u(x,y) as a function of x and y.
title('Approximate solution 9 pt Laplacian')


 x = gmres(A,b,100,10^-7);

% x = A\b;
% 
figure(3)
u = x;
ugrid = reshape(u,n-1,n-1);  % This makes the Nx1 vector u into an n-1 x n-1
                             % array.
surf([h:h:1-h]',[h:h:1-h]',ugrid)  %This plots u(x,y) as a function of x and y.
title('Exact solution')




for i = 1 : n
    
    for j = 1 : n
        
        G(i,j) = g(XM(i,j),YM(i,j));
    
    end
    
end

u = x;
ugrid = G ;%reshape(u,n-1,n-1);  % This makes the Nx1 vector u into an n-1 x n-1


% mesh(XM(1:40,1:40),YM(1:40,1:40),G)
% 
% surf([h:h:1-h]',[h:h:1-h]',ugrid)
% title('Exact solution')






function [Ax,G] = Lap_mult(x0)



g = inline('cos(x + 2*y)','x','y');


[n m] = size(x0);

n = sqrt(n);
h = 1/(n+1);

X = reshape(x0,[n,n]);

x = 0 : h : 1;
y = 0 : h : 1;

[XM,YM] = meshgrid(x,y);
G = [];
for i = 1 : n
    
    for j = 1 : n
        G(i,j) = g(XM(i,j),YM(i,j));
    
    end
    
end

omega = 1;

for i = 2 : n - 1
    for j = 2 : n - 1

   
        
    update = 1/6 * (X(i+1,j+1) + X(i+1,j-1) + X(i-1,j+1) + X(i-1,j-1)) + ...
            2/3 * (X(i+1,j) + X(i-1,j) + X(i,j+1) + X(i,j-1)) - 10/3 * X(i,j);
    X(i,j) = X(i,j) + 3/8 * omega * update;


    end 
end


n = n ^ 2;

Ax = reshape(X,[n 1]);

end


function [X,K,H,res_new,time,counter] = GMRES_1(x_0,A,b,m,tol)
 % to get CPU runtime (in miliseconds)
tic;
% initialize

 Ax = Lap_mult(x_0); 
 
 res = b - Ax;

beta = norm(res);
H = zeros(m - 1, m);
vec = res/beta;
n = length(A);
% call arnoldi orth function
[K,H,counter] = Arnoldi(vec,A,m,n,tol);

Y = (A * K(:,1:counter)) \ b
 %Y = Lap_mult(K(:,1:counter)) \ b;
% here are the solutions
X = x_0 + K(:,1:counter) * Y;
res_new = norm( A - b );
time = toc;
end

  
function [K,H,counter] = Arnoldi(vec,A,m,n,tol)

tol = 10 ^ -7;
% Initialize vectors K,H,W
 K = zeros(n,m);
%  H = zeros(m - 1,n);
counter = 0;
% Normalize the first
K(:,1) = vec / norm(vec);


    for j = 1:m
        
           w = Lap_mult(K(:,j));
         % w = A * K(:,j);
            
        for i = 1 : j 
            
            % the Hess1enberg matrix
            H(i,j) = K(:,i)' * w;
            w = w - H(i,j) * K(:,i);

        end
        
        
  
       % just to initialize the error
       if j < 2
           
            res = norm(w);
            
       end
       
       
        % here is if we get a lucky break
        
        if norm(w) == 0
            
                disp("Lucky break");
                return;
        end
        
        disp("Residual is :")
        disp(res)
        
        if res > tol
                            
                    H(j+1,j) = norm(w);
                
                
                y_m = inv(H'* H)* H';
                res = norm(w) * abs(y_m(j));
            
                K(:,j+1) = w / H(j+1,j);
        
        else
            
            return;
        
        end    
    
         [counter  bounter] = size(H);
          
    end
    
 
    
    
end
