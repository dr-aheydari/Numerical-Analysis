% Ali Heydari
% Iterative methods
% GMRES Method

n = 400;
A=full(gallery('tridiag',n,-1,2,-1));

b = zeros(n,1);
b(1) = 1;
% we will iterate up to this
m = n;
% arbitrary initial guess
x_0 = zeros(n,1);
n = size(x_0,1);

for i = 1:5
    tol = 10^-i;
[X,K,H,r_new,time,counter] = GMRES_1(x_0,A,b,m,tol);

fprintf("iter: %i    time:%i",counter,time);
end

function [X,K,H,res_new,time,counter] = GMRES_1(x_0,A,b,m,tol)
 % to get CPU runtime (in miliseconds)
tic;
% initialize
res = b - A * x_0;
beta = norm(res);
H = zeros(m - 1, m);
vec = res/beta;
n = length(A);
% call arnoldi orth function
[K,H,counter] = Arnoldi(vec,A,m,n,tol);

Y = (A * K(:,1:counter)) \ b;

% here are the solutions
X = x_0 + K(:,1:counter) * Y;
res_new = norm( A * X - b );
time = toc;
end

  
function [K,H,counter] = Arnoldi(vec,A,m,n,tol)

tol = 10 ^ -5;
% Initialize vectors K,H,W
 K = zeros(n,m);
%  H = zeros(m - 1,n);
counter = 0;
% Normalize the first
K(:,1) = vec / norm(vec);


    for j = 1:m
       
         w = A * K(:,j);
            
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
        
        if res > tol
                            
                    H(j+1,j) = norm(w);
                
                
                y_m = inv(H'* H)* H';
                res = norm(w) * abs(y_m(j));
            
                K(:,j+1) = w / H(j+1,j);
        
        else
            
            return;
        
        end    
    
         [counter  bounter] = size(H);
%         disp(bounter);

    end
end

