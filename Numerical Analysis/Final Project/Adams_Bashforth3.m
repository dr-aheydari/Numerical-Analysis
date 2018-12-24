    % Ali Heydari
    % LMM 2(b)
    % Math 231 Proj 6

    

    % to get CPU time

    % desired interval

    a = 2;
    b = 3;

    % size of each step
    % h_step = [0.1 0.05];


    h_step = [0.1 0.01 0.025 0.05 0.005 0.001 ] ;

    % initial value (we have w_1 and w_2 from the Euler's method approximation


    w_0 = [1.0000];


    % the function being solved

    f{1} = @(t,y) 1 + (t - y)^2;
    f{2} = @(t,x) -t * x;
    % f{3} = @(t,x) exp(-2 * t) - 2 * x;


    % To make our life a lot easier with the plots

    str = ["y' = 1 + (t - y)^2 "];


    % function call for h = 0.1

    [w, error_last,cpu_time] = adam_bash_3(a,b,h_step,w_0,f,str);

    % % function call for h = 0.05



    %
    % ylabel(" Log of CPU Time");
    % title("Precision Diagram for Adam-Bashforth 3rd order vs Euler's");











    function [w,error_last,cpu_time] = adam_bash_3(a,b,h_step,w_0,f,str)

    % Number of iterations


    z = 1;


    % v will be the cell storing the approximated solutions

    for k = 1 : length(h_step);

        h = h_step(k);

        N = (b - a) / h;

        % initialization
        w = zeros(N,1);
        t = zeros(N,1);
        v = cell(3,1);


        for i = 1 : 3

            v{i} = zeros (N,1);

        end





        for i = 1 : 1

            % for easier use
            g = f{i};
            w = v{i};
            % initial value, given to us
            w(1) = w_0(i);
            u(1) = w_0(i);
            q(1) = w_0(i);
            t(1) = 2;
            t(2) = 2.01;
            t(3) = 2.02;
            d = zeros(1,N);
            u = v{2};
            q = v{3};

            % ODE 45
            %             N = (b - a) / h
            t = 2 : h : 3 ;
            [time sol] = ode45(g,t,w_0);

            for m = 2:4

                K1 = h * g(t(m-1),w(m-1));
                K2 = h * g(t(m-1) + h/2, w(m-1) + K1 /2);
                K3 = h * g(t(m-1) + h/2, w(m-1) + K2 /2);
                K4 = h * g(t(m-1) +h, w(m-1) + K3);

                w(m) = w(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                u(m) = u(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;
                q(m) = q(m-1) + (K1 + 2 * K2 + 2 * K3 + K4) * 1/6;

            end

            tic;
            for j = 3 : N

                if i == 1

                    t(j + 1) = a + h * j;

                end




                % Adams-Bashforth three step
                w(j + 1) = w(j) + h/12 * (23 * g(t(j),w(j)) - 16 * ...
                    g(t(j-1),w(j-1)) + 5 * g(t(j-2),w(j-2)));

%                 sol(j)
%                 w(j)
                error_last(z,1) = norm(sol(j) - w(j));

            end
            cpu_time(z,1) = toc;
            % BDF 3 step

            tic;
            for j = 3 : N


                u(1) = w(1);
                u(2) = w(2);
                u(3) = w(3);


                % Backwards differenciation Formula-Three step
                t_s = @(i) a + (i-1)*h;
                fp = @(uFixPoint,u,g,c) (1/11) * (18 * u(j) - 9 * u(j-1) + 2 * u(j-2) ...
                    + 6 * h * g(t_s(c + 1), uFixPoint) );

                uFixPoint = u(j);

                while (abs(uFixPoint - fp(uFixPoint,u,g,j))> 1e-8)

                    uFixPoint = fp(uFixPoint,u,g,j);


                end

                u(j+1) = uFixPoint;

                error_last(z,2) = norm(sol(j) - u(j));


                %              u(j+1) = (1/11) * (18 * u(j) - 9 * u(j-1) + 2 * u(j-2)  ...
                %                  + 6 * h * w(j+1));


            end

            cpu_time(z,2) = toc;

            % Adams Moulton 2-step method:


            tic;
            for j = 3 : N


                % Adams Moulton 2-step method:

                q(1) = w(1);
                q(2) = w(2);
                q(3) = w(3);


                fp = @(wFixPoint,q,g,c) q(c) + (h*(1/12)) * (5 * g(t_s(c + 1), wFixPoint)...
                    + 8 * g(t_s(c), q(c)) - g(t_s(c - 1), w(c - 1)));

                wFixPoint = q(j);

                while (abs(wFixPoint - fp(wFixPoint,q,g,j))> 1e-8)

                    wFixPoint = fp(wFixPoint,q,g,j);


                end

                q(j+1) = wFixPoint;

                error_last(z,3) = norm(sol(j) - q(j));


            end

            cpu_time(z,3) = toc;

            % storing the current solution
            v{1} = w;
            v{2} = u;
            v{3} = q;


            z = z + 1;

        end

        % initialization





        figure(k)
%         size(t)
%         disp("SOL")
%         size(w)
        hold on
        plot(time,sol,'b-o');
        plot (t,w,'r-*');
        plot (t,u,'g-x');
        plot (t,q,'k-+');

        disp_h = sprintf('Step size (h) = %f',h);
        
    xlabel (" Time ");
    ylabel (" y(t) ");
    title ({'Approximate Solution vs. Actual Solution for : ';  str(i); ...
        disp_h});
    legend("Actual Solution", "Adams-Bashforth 3 step",...
        "BDF 3 step", "Adams-Moulton 2 step");



    end


    hold off


    % returning the final approximations
    w = v;

    end


