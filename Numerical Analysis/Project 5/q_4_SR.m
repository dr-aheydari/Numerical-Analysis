% Ali Heydari
% Math 231 Homework 4
% SR Method (Q4)

%% Interactive Interface (with user input)
% get input
% A = input(?Please enter the matrix A: ?);
% x_init = input(' Please enter a vector of initial guesses for x :');
% delta = input(?Please enter the desired tolerance: ?);

%% Non Interactive (without user input)

% making the tridiagonal matrix
n = 10;
A=full(gallery('tridiag',n,-1,2,-1));

b = zeros(n,1);
b(1) = 1;

% arbitrary initial guess
x_init = zeros(n,1);
n = size(x_init,1);
error = 1;
w = 1.5;

% different norms

% l2 norm
L = 2;

%% Different norms and Tolerances
% find the solution based on what norm and what error:
% delta = 10^-2;
% [x_sol_l2_1 ,counter_l2_1] = gauss_seidel(A, x_init, delta, error, n, b,...
%     L, w);

delta = 10^-4;
[x_sol_l2_2 ,counter_l2_2] = gauss_seidel(A, x_init, delta, error, n, b,...
 L, w);

% l_infty norm (denoted as l8)
% L = Inf;
% delta = 10^-2;
% [x_sol_l8_1 ,counter_l8_1] = gauss_seidel(A, x_init, delta, error, n, b,...
%     L, w);
% 
% delta = 10^-5;
% [x_sol_l8_2 ,counter_l8_2] = gauss_seidel(A, x_init, delta, error, n, b,...
%     L, w);


%% Output Formatting (Credit to Shayna for showing me this)

% Col1 = Nr1(1:end-1);
% Col2=abs(Nr1(2:end)-Nr1(1:end-1));
% Col3=abs(Col1-pN1);
% T=table(Col1',Col2',Col3')


%% Function 

function[x_sol,counter] = gauss_seidel(A, x_init, delta, error, n, b, L, w)
counter = 0;

% the algorithm
while error > delta
    
    x_old = x_init;
    
    for i = 1:n
        
        % this is essentially a place holder for the sum part
        AX = 0;
        
        % since i not= j, we do two for loops
        for j = 1 : i-1
            
            AX =  (x_old(j) * A(i,j)) + AX ;
            
        end
        
        for j = i+1 : n
            
            AX = (x_old(j) * A(i,j)) + AX ;
            
        end
        
        % now get the value of each
        x_init(i) = (1 / A(i,i)) * (b(i) - AX);
        
    end
    
    % update error and counter 
    
    if L == Inf
    
        error = norm((x_old - x_init),Inf);
    
    end
    
    if L == 2
    
            error = norm((x_old - x_init));        
    end
    
    counter = counter + 1;
end

% storing our answers 
counter = counter - 1;
x_sol = x_init;

end

%% Output