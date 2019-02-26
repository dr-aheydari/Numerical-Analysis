% Ali Heydari
% Getting a reference solution
% Math 231: Final Project

k = 2:.01:3;

f = @(t,y) 1 + (t - y)^2
y_2 = 1;

% to set the absolute and relative tolerance to 10^-13
odeset('RelTol', 1e-13, 'AbsTol', 1e-13);


[time sol] = ode45(f,k,y_2);


plot(time,sol,'r-o');
xlabel("Time");
ylabel("Approximated y(t)")
title("Approximated solutions using multi-step methods");