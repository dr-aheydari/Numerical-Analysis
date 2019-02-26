% Ali Heydari
% Math 231, hw3
% Bisection Method


%%% Interactive Interface (with user input)
% % get input 
% a = input('Please enter a value for the lower bound a: ');
% b = input('Please enter a value for the upper bound (b) : ');
% delta = input('Please enter the desired tolerance: ');
% f = input('Please enter f(x)?(type @(x) [then the function] ');





% see if any of the boundaries are a root 

%%% Non Interactive (without user input)

f = @(x) x^3 - 3*x + 2;


retur = 0;
counter = 0;
x_k = ones(1,10);
error = zeros(1,10);
e_n = zeros(1,10);

delta = 10^-6
a = -4;
b = 0;


%%% Method
fa = f(a);

if fa == 0 
    
    root = a; 
    retur = 1;
    
end; 

fb = f(b);

if fb == 0 
    root = b; 
    retur = 1; 

end; 

% if the boundaries are not the root then do bisection

if retur ~= 1

    % check if the user hasnt lost their mind
    if sign(fa) == sign(fb)
       
        display('Error: f(a) and f(b) have same sign.')
        
        retur = 1;
       
    end;
         
    % if all is gucci
    
    
           while abs(b-a) > delta && retur ~= 1
                % As Mayya said in class, keep going until Iterate I <= 2delta
                
                counter = counter + 1;
                
                c = (a+b)/2;
                fc = f(c);
                
                if fc == 0 
                        root = c; 
                        retur = 1;
                end;
                
                % cehk to see which side of the interval we want
                
                if sign(fc) == sign(fa)
                        a = c;
                        fa = fc;
                else 
                        b = c;
                        fb = fc;
                end;
                
            end;
            
            % hopefully we got want we need
            root = (a+b)/2;
end

fprintf('The found root is: %i \n',root);
fprintf('Total iterations: %i \n', counter);


%%% Output Formatting

disp(" ");
disp(" ");
fprintf('The root of the function is at x = %i \n', root);
fprintf('Number of iterations: %i \n', counter);
disp(" ");
disp(" ");

disp("        pn              |p_{n+1} - p_n|       e_n = |pn - p| ")
for i= 1 : counter
    
fprintf("%i   %i         %i         %i\n",i ,x_k(i),error(i),e_n(i));

end

%%% Outputs

