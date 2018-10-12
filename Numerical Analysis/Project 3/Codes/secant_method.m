% Ali Heydari
% Math 231, hw3
% Secant Method


% Problem 3)part a) second root x = 1

x_k = zeros(1,10);
error = zeros(1,10);
e_n = zeros(1,10);

%%% Interactive Interface (with user input)

% % get initial conditions
% x_k(1) = input('Please enter the initial guess x_0: ');
% x_k(2) = input('Please enter the initial guess x_1: ');
% delta = input('Please enter the desired tolerance: ');
% f = input('Please enter f(x)?(type @(x) [then the function] ');

%%% Non Interactive (without user input)
f = @(x) x^3 - 3*x + 2;
% our two initial guesses
x_k(1) = 1.2;
x_k(2) = 0;
xk = x_k(1); 
fx = f(x_k(1));
delta = 10^-6 ;

% evaluate function at x_o: f(x_o)
fx = f(x_k(2)); 

counter = 2; % counter

% keep finding the root until f(x_k) = 0
while abs(f(x_k(counter))) > delta
    
% secant method formula
x_k(counter+1) = x_k(counter) - (f(x_k(counter))*(x_k(counter) ...
    - x_k(counter-1)))/( f(x_k(counter)) - f(x_k(counter-1)));
% %display xk and f(xk)
 xk = x_k(counter+1);
 fx = f(x_k(counter+1));
% find error error:

e_n(counter) = abs(xk - (-2));
     if fx <= 1 * 10 ^ -6
            
            root = xk;   
     end

error(counter) = abs(x_k(counter+1) - x_k(counter));
counter = counter + 1;

end

%%% Output Formatting

disp(" ");
disp(" ");
fprintf('The root of the function is at x = %i \n', root);
fprintf('Number of iterations: %i \n', counter);
disp(" ");
disp(" ");

disp("        pn              |p_{n+1} - p_n|       e_n = |pn - p| ")
for i= 1 : counter - 1
    % to only output 10 values
    if i < 11
        fprintf("%i   %i         %i         %i\n",i ,x_k(i),error(i),e_n(i));
    end
end

%%% Outputs