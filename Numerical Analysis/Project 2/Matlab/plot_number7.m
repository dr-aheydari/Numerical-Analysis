clear('all');
x = -3 : 0.5 : 5;

hold on

y = x.^2/2 + 5*x/2 + 1;

p1 = plot(x,y);
set(p1,'LineWidth',4);
xlabel('x');
ylabel('y');

p2 = plot(0,1,'O');
set(p2,'LineWidth',4);

t1 = text(0,2,'\downarrow (0,1)');
set(t1,'LineWidth',4);
p3 = plot(1,4,'O');
set(p3,'LineWidth',4);

text(1,5,'\downarrow (1,4)')
p4 = plot(2,8,'O');
set(p4,'LineWidth',4);

text(2,9,'\downarrow (2,8)')
legend('f(x) = x^2/2 + 5x/2 + 1','(0,1)','(1,4)','(2,8)');
hold off