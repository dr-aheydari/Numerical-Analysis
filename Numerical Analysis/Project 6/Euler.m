% Ali Heydari
% Euler's Method
% Math 231 Proj 6

%% non-interactive

% desired interval

a = 0;
b = 10;

% size of each step
h = 0.1;

% initial value

w_0 = [1 1 1/10];

% the function being solved

f{1} = @(t,x) t^2 - x;
f{2} = @(t,x) -t * x;
f{3} = @(t,x) exp(-2 * t) - 2 * x;

% given solutions
y = cell(3,1);

y{1} = @(t) -1 * exp(-t) + t^2 - 2 * t + 2;
y{2} = @(t) exp(-t^2 / 2);
y{3} = @(t) 1/10 * exp(-2 * t) + t * exp(-2 * t);

% To make our life a lot easier with the plots

str = ["y' = t^2 - y ","y' = -t * y","y' = e^{-2t} - 2y"];
      

% function call

[w,error] = euler(a,b,h,w_0,f,y,str);
    


function [w,error] = euler(a,b,h,w_0,f,y,str)

% Number of iterations

N = (b - a) / h;

w = zeros(N,1);

t = zeros(N,1);

v = cell(3,1);

for i = 1 : 3

    v{i} = zeros (N,1);
    
end


for i = 1 : 3

    g = f{i};
    w = v{i};
    w(1) = w_0(i);
    
    for j = 1 : N
    
        w(j + 1) = w(j) + h * g(t(j),w(j));
        size(w)
        
        if i == 1
            
             t(j + 1) = a + h * j ;
        
        end
        
    end
    
        v{i} = w;
   
    
end

error = zeros(N,3);
    A = zeros(N,1);

for  i = 1 : 3
    
    g = y{i};
    w = v{i};

    for j = 1 : N
        
        A(j) = g(t(j));
        
        error(j,i) = w(j) - g(j);
        
    end
    
    figure(i);
    hold on
    plot (t(2:101), A ,'b--o');
    plot (t(2:101),w(1:100),'r');
    
    xlabel (" Time ");
    ylabel (" y(t) ");
    title ({'Approximate Solution vs. Actual Solution'; str(i)});
    legend("Actual Solution", "Approximation");
   
    
end

 hold off

w = v;

end