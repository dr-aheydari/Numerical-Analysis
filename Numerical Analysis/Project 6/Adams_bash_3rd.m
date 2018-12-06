% Ali Heydari
% Adams-Bashforth Method
% Math 231 Proj 6

%% Adam's Bashforth

% to get CPU time

% desired interval

a = 0;
b = 10;

% size of each step
% h_step = [0.1 0.05];

for ii = 1 : 20

h_step = 0.5 / ii
h(ii) = h_step;

boo = tic;

% initial value (we have w_1 and w_2 from the Euler's method approximation

w_0 = [1 1 1/10];
w_1 = [0.900000000000000 1  0.180000000000000];
w_2 = [0.811000000000000 0.990000000000000 0.225873075307798];

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
      

% function call for h = 0.1

[w] = adam_bash_3(a,b,h_step,w_0,w_1,w_2,f,y,str);
    
% % function call for h = 0.05
% h = 0.05
% [w,error] = euler(a,b,h,w_0,f,y,str);
     
e(ii) = toc(boo);

end

% 
hold on
p = scatter(log(h),log(e),'X');

q = polyfit(h,e,2)
y = 0 + q(1,1)*h.^2 + q(1,2)*h + q(1,3); %q(1,1)*h.^3
plot(log(h),log(y));

% ylabel(" Log of CPU Time");
% title("Precision Diagram for Adam-Bashforth 3rd order vs Euler's");

%% Euler's Method

% desired interval



%  clear("all");

for k = 1 : 20
a = 0;
b = 10;
% size of each step

h_step = 0.5 / k
h(k) = h_step;
aoo = tic;
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
      

% function call for h = 0.1

[w] = euler(a,b,h_step,w_0,f,y,str);
    
% % function call for h = 0.05
% h = 0.05
% [w,error] = euler(a,b,h,w_0,f,y,str);

l(k) = toc(aoo);
end

% hold on
% p = scatter(h,e,'X');
% 
% q = polyfit(h,e,3)
% y = q(1,1)*h.^3 + q(1,2)*h.^2 + q(1,3)*h + q(1,4);
% plot(h,y);
% 
% 
% hold on
% figure(1)
% p = scatter(log(h),log(e),'X');
% q = polyfit(h,e,2)
% y = 0 + q(1,1)*h.^2 + q(1,2)*h + q(1,3); %q(1,1)*h.^3
% plot(log(h),log(y));


j = scatter(log(h),log(l),'o');
z = polyfit(h,l,2)
v = z(1,1)*h.^2 + z(1,2)*h.^1 + z(1,3)*h.^0; % + z(1,4);
plot(log(h),log(v));

xlabel("Log of h (step size)");
ylabel("Log of CPU Time");
title("Precision Diagram for Adam-Bashforth 3rd order vs Euler's");
legend("Adams-Bashforth Method", "Polynomial Fit for AB", ...
    "Euler's Method", "Polynomail fit for Euler's");

function [w] = euler(a,b,h_step,w_0,f,y,str)

% for k = 1 : 2
% h = h_step(k);

h = h_step;
% Number of iterations


N = (b - a) / h;

% initialization
w = zeros(N,1);
t = zeros(N,1);
v = cell(3,1);

% c = 17840;
% 
% for i = 1 : 10000
%     
%     
%     c = 1000 * 1000 * c;
% end
% v will be the cell storing the approximated solutions
for i = 1 : 3

    v{i} = zeros (N,1);
    
end

h = h_step;

for i = 1 : 3

    % for easier use 
    g = f{i};
    w = v{i};
    % initial value, given to us
    w(1) = w_0(i);
    
    for j = 1 : N
    
        % Euler's 
        w(j + 1) = w(j) + h * g(t(j),w(j));
        
        % to save computaion time, since all the t's will be the same
        % for the consecutive i's
        if i == 1
            
             t(j + 1) = a + h * j ;
        
        end
        
    end
    
    % storing the current solution
        v{i} = w;
   
    
end

% initialization
error = zeros(N,3);
    A = zeros(N,1);

for  i = 1 : 3
    
    % for easier use 
    g = y{i};
    w = v{i};

    for j = 1 : N
        
        % evaluating the exact solution at each time
        A(j) = g(t(j));
        
        % computing the error
        error(j,i) = norm(w(j) - A(j));
       
        % to get the error at time t == 2
        
%         if h == 0.1
%        
%             if j == 20
%            
%                 fprintf("For the differential equation %s \n", str(i));
%                 fprintf("The error at time = 2 is : %i for h = %i \n",...
%                     error(j,i),h);
%                 disp(" ");
%             end
%         end
%         
%          if h == 0.05
%                 
%                 if j == 40
%                  
%                 fprintf("For the differential equation %s \n", str(i));
%                 fprintf("The error at time = 2 is : %i for h = %i \n",...
%                     error(j,i),h);
%                 disp(" ");
%         
%             
%                 end
%                 
%          end
%         
%     end
    
%    plotting stuff
%     figure(i);
%     hold on
%     plot (t(2:101), A ,'b--o');
%     plot (t(2:101),w(1:100),'r');
  
%     if k == 1
%         
%         plot (t(2:101), A ,'b--o');
%         plot (t(2:101),w(1:100),'r');
%     
%     end
%     
%     if k == 2
%       
%        plot (t(1:200),w(1:200),'g');    
%     
%     end
%     
%     
%     
%     
%     xlabel (" Time ");
%     ylabel (" y(t) ");
%     title ({'Approximate Solution vs. Actual Solution for : ';  str(i)});
%     legend("Actual Solution", "Approximation for h = 0.1",... 
%         "Approximation for h = 0.05");
%    
    
   
    
%     % the error plots
%       figure(3 + i)
%      
%       hold on
%       title ({'Plot of the Local Truncation Error vs h for:';  str(i)})
%      error_full{i,k} = error(:,i);
%      err{i,k} = norm(error(:,i));
% %      plot(h_step(k),log(norm(error(:,i))),'-x');
%      
%      er(k) = err{i,k};
%      plot(h_step(k),er(k),'b--o');
%      
%       if k == 2
%          
%          p{i} = polyfit(h_step,er,1);
%          w = p{i};
%           x = 0:0.1:.2;
%           if i == 1
%               % i had a small bug, so I mannually calculated the line for 
%               % this one
%               
%               c = .8845 + 5.186*(x - 0.1);
%               
%           else
%               
%               c = w(1)*x + w(2);
%           end
%          
%           plot(x,c)
%          
%       end
% 
%        xlabel (" h (Step Size) ");
%        ylabel (" Local Truncation Error ");
%       
%       
%       hold off
end



% end


% returning the final approximations
w = v;

end


end






function [w] = adam_bash_3(a,b,h_step,w_0,w_1,w_2,f,y,str)

% Number of iterations



% v will be the cell storing the approximated solutions

 for k = 1:2
    
    h = h_step;

    N = (b - a) / h;

% initialization
w = zeros(N,1);
t = zeros(N,1);
v = cell(3,1);


for i = 1 : 3
    
    v{i} = zeros (N,1);
    
end


for i = 1 : 3

    % for easier use 
    g = f{i};
    w = v{i};
    % initial value, given to us
    w(1) = w_0(i);
    w(2) = w_1(i);
    w(3) = w_2(i);
    
    for j = 3 : N
    
        % Euler's 
        w(j + 1) = w(j) + h/12 * (23*g(t(j),w(j)) - 16*g(t(j-1),w(j-1))+ ...
            5*g(t(j-2),w(j-2)));
        
        % to save computaion time, since all the t's will be the same
        % for the consecutive i's
        if i == 1
            
             t(j + 1) = a + h * j ;
        
        end
        
    end
    
    % storing the current solution
        v{i} = w;
   
    
end

% initialization
error = zeros(N,3);
    A = zeros(N,1);

for  i = 1 : 3
    
    % for easier use 
    g = y{i};
    w = v{i};

    for j = 1 : N
        
        % evaluating the exact solution at each time
        A(j) = g(t(j));
        
        % computing the error
        error(j,i) = norm(w(j) - A(j));
        
        % to get the error at time t == 2
        
%         if h == 0.1
%        
%             if j == 20
%            
%                 fprintf("For the differential equation %s \n", str(i));
%                 fprintf("The error at time = 2 is : %i \n",error(j,i));
%                 disp(" ");
%             end
%         end
%         
%          if h == 0.05
%                 
%                 if j == 40
%                  
%                 fprintf("The error at time = 2 is : %i \n ",error(j,i));
%         
%             
%                 end
%                 
%          end
%         
%     end
%     
    
    % plotting stuff
        figure(i);
        hold on
    
    
    
    if k == 1
        
        plot (t(2:101), A ,'b--o');
        plot (t(2:101),w(1:100),'r');
    
    end
    
    if k == 2
      
       plot (t(1:200),w(1:200),'g');    
    
    end
    
    
    xlabel (" Time ");
    ylabel (" y(t) ");
    title ({'Approximate Solution vs. Actual Solution for : ';  str(i)});
    legend("Actual Solution", "Approximation for h = 0.1",... 
        "Approximation for h = 0.05");
   
    the error plotsclear
    figure(3 + i)
    loglog( error(:,i));
    title ({'Log-Log Plot of the Error for:';  str(i)})
    
    
end


%  hold off

end
% returning the final approximations
w = v;

end