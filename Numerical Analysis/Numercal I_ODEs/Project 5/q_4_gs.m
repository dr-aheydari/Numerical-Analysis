% Ali Heydari
% Math 231 Homework 4
% Gauss_Seidel Method (Q4) N = 200

%% Interactive Interface (with user input)
% get input
% A = input(?Please enter the matrix A: ?);
% x_init = input(' Please enter a vector of initial guesses for x :');
% delta = input(?Please enter the desired tolerance: ?);

%% Non Interactive (without user input)

% making the tridiagonal matrix
n = 400;
A=full(gallery('tridiag',n,-1,2,-1));

b = zeros(n,1);
b(1) = 1;

% arbitrary initial guess
x_init = zeros(n,1);
n = size(x_init,1);
error = 1;


% different norms

% l2 norm
L = 2;

%% Different norms and Tolerances
% find the solution based on what norm and what error:
% delta = 10^-2;
% [x_sol_l2_1 ,counter_l2_1] = gauss_seidel(A, x_init, delta, error, n, b, L);
for i = 1 : 5
delta = 10^-i;

[x_sol_l2_2 ,counter_l2_2,time] = gauss_seidel(A, x_init, delta, error, n, b, L);

fprintf("iter: %i    time:%i\n",counter_l2_2,time);
end
% l_infty norm (denoted as l8)
% L = Inf;
% delta = 10^-2;
% [x_sol_l8_1 ,counter_l8_1] = gauss_seidel(A, x_init, delta, error, n, b, L);
% 
% delta = 10^-5;
% [x_sol_l8_2 ,counter_l8_2] = gauss_seidel(A, x_init, delta, error, n, b, L);
% 

%% Output Formatting 

% fprintf(" APPROXIMATIONS: \n")
% disp(" ")
% fprintf(" For l-2 norm and 10^-2   |   For l-2 norm and 10^-5 | \n ");
% disp("----------------------------------------------------- ")
% for i = 1 : size(x_sol_l2_1)
%     
%     fprintf("        %f          |         %f\n", x_sol_l2_1(i),x_sol_l2_2(i));
% 
% end
% disp("----------------------------------------------------- ")
% fprintf(" NUMBER OF ITERATIONS\n ")
% disp(" ")
% fprintf("For l-2 norm and 10^-2 | For l-2 norm and 10^-5 |\n");
% disp("----------------------------------------------------- ")
%  fprintf("        %f          |         %f\n", counter_l2_1,counter_l2_2);
% disp(" ")
% fprintf(" APPROXIMATIONS: \n")
% fprintf(" For l-Inf norm and 10^-2  |  For l-Inf norm and 10^-5 | \n ");
% disp("----------------------------------------------------- ")
% for i = 1 : size(x_sol_l8_1)
%     
%     fprintf("        %f           |         %f\n", x_sol_l8_1(i),x_sol_l8_2(i));
% 
% end
% 
% disp("----------------------------------------------------- ")
% disp(" ")
% fprintf(" NUMBER OF ITERATIONS\n ")
% fprintf("For l-Inf norm and 10^-2 | For l-Inf norm and 10^-5 |\n");
% disp("----------------------------------------------------- ")
%  fprintf("        %f          |         %f\n", counter_l8_1,counter_l8_2);


%% Function 

function[x_sol,counter,time] = gauss_seidel(A, x_init, delta, error, n, b, L)
tic;

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
time = toc;
end

%% Output