    % Ali Heydari
    % LMM, 2(c)
    % Math 231 Proj 6

    %% Adam's Bashforth

    % to get CPU time

    % desired interval

    a = 0;
    b = 10^-3;

    % size of each step
    % h_step = [0.1 0.05];


    h_step = [0.01]; %0.01 0.025 0.05 0.005 0.001 ] ;

    % initial value (we have w_1 and w_2 from the Euler's method approximation


    w_0 = [2.0000 0];


    % the function being solved

    f{1} = @(y) y ;
    f{2} = @(x,y) 500 * (1 - x^2) * y - x; 
    
    


    % To make our life a lot easier with the plots

    str = ["y'_1 = y, y'_2 = 500 * (1 - y_1 ^ 2) * y_2 - y_1"];


    % function call for h = 0.1

    [w, error_last,cpu_time] = adam_bash_3(a,b,h_step,w_0,f,str);

    % % function call for h = 0.05



    %
    % ylabel(" Log of CPU Time");
    % title("Precision Diagram for Adam-Bashforth 3rd order vs Euler's");











    function [w,sol,cpu_time] = adam_bash_3(a,b,h_step,w_0,f,str)

    % Number of iterations


    z = 1;


    % v will be the cell storing the approximated solutions

    for k = 1 : length(h_step);

        h =  0.0250 * 10^-3 

        N = (b - a) / h

        % initialization
        w = zeros(N + 1,1);
        t = zeros(N,1);
        v = cell(3,1);


        for i = 1 : 3

            v{i} = zeros (N,1);

        end





        for i = 1 : 1

            % for easier use
            g1 = f{1}
            g2 = f{2}
            
            w1 = v{1};
            w2 = v{2};
            % initial value, given to us
            w1(1) = w_0(1);
            w2(1) = w_0(2);
            w1(2) = w_0(1);
            w2(2) = w_0(2);
            w1(3) = w_0(1);
            w2(3) = w_0(2);
 
            
             u1 = v{2};
             u2 = v{2};
             q1 = v{3};
             q2 = v{3};
             
             u1(1) = w_0(1);
             u2(1) = w_0(2);
             q1(1) = w_0(1);
             q2(1) = w_0(2);

            % ODE 45
            %             N = (b - a) / h
             t = 0 : h : 10^-3 ;
            
           [time sol] = ode45(@vdp1,[0 10^-3],[2; 0])
%             [time sol] = ode45(g,t,w_0);
% 
            for m = 2:4

                K1 = h * g1(w1(m-1));
                K2 = h * g1( w1(m-1) + K1 /2);
                K3 = h * g1( w1(m-1) + K2 /2);
                K4 = h * g1( w1(m-1) + K3);

                w1(m) = w1(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                
                
                
                
                
                 K1 = h * g2(t(m-1),w2(m-1));
                 K2 = h * g2(t(m-1) + h/2, w2(m-1) + K1 /2);
                 K3 = h * g2(t(m-1) + h/2, w2(m-1) + K2 /2);
                 K4 = h * g2(t(m-1) +h, w2(m-1) + K3);
% 
                 w2(m) = w2(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                 u1(m) = u1(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                 q1(m) = q1(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                
                 u1(m) = u1(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                 q1(m) = q1(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;

            end





            tic;
            for j = 3 : N 

                if i == 1

                    t(j + 1) = a + h * j;

                end




                % Adams-Bashforth three step
                  w1(j + 1) = w1(j) + h/12 * (23 * g1(w2(j)) - 16 * ...
                  g1(w2(j-1)) + 5 * g1(w2(j-2)));
                
%              
                  w2(j + 1) = w2(j) + h/12 * (23 * g2(w1(j),w2(j)) - 16 * ...
                      g2(w1(j-1),w2(j - 1)) + 5 * g2(w1(j-2),w2(j - 2)));

%                 sol(j)
                 

                 error_last(z,1) = norm(sol(j) - w(j));

            end
            
            w(:,1) = w1(:);
            w(:,2) = w2(:);
            
            
            cpu_time(z,1) = toc;
            
            
            % BDF 3 step

            tic;
            for j = 3 : N


                % Backwards differenciation Formula-Three step
                t_s = @(i) a + (i-1)*h;
%                 fp = @(uFixPoint,u,g1,g2,c) (1/11) * (18 * u(j) - 9 * u(j-1) ...
%                     + 2 * u(j-2) + 6 * h * g2(uFixPoint,uFixPoint));


                fp1 = @(uFixPoint,u1,g1,c) (1/11) * (18 * u1(j) - 9 * u1(j-1) ...
                     + 2 * u1(j-2) + 6 * h * g1(uFixPoint));
                 
                


                u1FixPoint = u1(j);
                u2FixPoint = u2(j);

                while (abs(u1FixPoint - fp1(u1FixPoint,u1,g1,j)) > 10^-2)

                    u1FixPoint = fp1(u1FixPoint,u1,g1,j);


                end

                u1(j+1) = u1FixPoint;
                
                fp2 = @(uFixPoint,u2,g2,c) (1/11) * (18 * u2(j) - 9 * u2(j-1) ...
                     + 2 * u2(j-2) + 6 * h * g2(u1(j+1),uFixPoint));
                
                
                while (abs(u2FixPoint - fp2(u2FixPoint,u2,g2,j)) > 10^-2)

                    u2FixPoint = fp2(u2FixPoint,u2,g2,j);


                end
                
                u2(j+1) = u2FixPoint;


% 
                 u = [];
                 u(:,1) = u1(:);
                 u(:,2) = u2(:);
%                 error_last(z,2) = norm(sol(j) - u(j));

% 
%                 u1(j+1) = (1/11) * (18 * u1(j) - 9 * u1(j-1) + 2 * u1(j-2)  ...
%                        + 6 * h * w(j+1));
%                    
%                     u(j+1) = (1/11) * (18 * u(j) - 9 * u(j-1) + 2 * u(j-2)  ...
%                        + 6 * h * w(j+1));


            end
% 
%             cpu_time(z,2) = toc;
% 
%             % Adams Moulton 2-step method:
% 
% 
%             tic;
 counter = 1;
            for j = 3 : N


                % Adams Moulton 2-step method:
                t_s = @(i) a + (i-1)*h;
                fp = @(qFixPoint,q,g1,c) q(c) + (h * (1/12)) * (5 * g1(qFixPoint)...
                    + 8 * g1(q(c)) - g1(q(c - 1)));
                
                q1FixPoint = q1(j);
                q2FixPoint = q2(j);
                
               

                while (abs(q1FixPoint - fp(q1FixPoint,q1,g1,j)) > 10^-2)

                    q1FixPoint = fp(q1FixPoint,u1,g1,j);

%                 counter = counter + 1 

                end

                q1(j+1) = q1FixPoint;
                
                fp = @(qFixPoint,q1,q2,g2,c) q2(c) + (h * (1/12)) * (5 * ...
                    g2(q1(j+1),qFixPoint)...
                    + 8 * g2(q1(c),q2(c)) - g2(q1(c-1),q2(c - 1)));
                
                while (abs(q2FixPoint - fp(q2FixPoint,q1,q2,g2,j)) > 10^-8)

                    q2FixPoint = fp(q2FixPoint,q1,q2,g2,j);


                end

                q2(j+1) = u2FixPoint;

                q = [];
                q(:,1) = q1(:);
                q(:,2) = q2(:);

                error_last(z,3) = norm(sol(j) - q(j));


            end

            cpu_time(z,3) = toc;
% 
%             % storing the current solution
%             v{1} = w;
%             v{2} = u;
%             v{3} = q;
% 
% 
%             z = z + 1;
% 
%         end
% 
%         % initialization
% 
% 
% 
% 
% 
%         figure(k)
% %         size(t)
% %         disp("SOL")
% %         size(w)
          hold on
           plot(time,sol(:,1),'b-o');
           plot(time,sol(:,2),'k-o');
           disp(size(t));
           plot (t,q(:,1),'r-*');
           plot (t,q(:,2),'g-*');
%          plot (t,q,'k-+');
% 
         disp_h = sprintf('Step size (h) = %f',h);
%         
    xlabel (" Time ");
    ylabel (" y(t) ");
    title ({'Approximate Solution vs. Actual Solution for : ';  str(i); ...
        disp_h});
   legend( {"Actual Solution: y_1","Actual Solution: y_2" ...
        "Adam-Moulton 2 Step: y_1", "Adam-Moulton 2 Step: y_2" ...
        },'FontSize',22);

% "BDF 3 step", "Adams-Moulton 2 step"

    end


%     hold off

    end
    % returning the final approximations
%     w = v;

    end

    function dydt = vdp1(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

dydt = [y(2); 500 * (1-y(1)^2)*y(2)-y(1)];
    end

