% Ali Heydari
% Graphing Region of Absolute Convergance
% Math 231 Final Proj 
% Reference: 
% https://math.boisestate.edu/~wright/courses/m566/StabilityDomains/html/StabilityDomains.html

LW = 'LineWidth'; lw = 1;
N = 1000;
th = linspace(0,2*pi,N);
r = exp(1i*th);

% this part is changed based on the specific methid
% f = @(r) 11*(r.^3-(18/11)*r.^2 + (9/11)* r - (2/11))./(6* r.^3);
f = @(r) (r.^2 - r) ./ r.^1
xi = f(r);
plot(xi,'k-',LW,lw),
hold on
fill(real(xi),imag(xi),'c')
x = -4.5: .2 : 7.5;
y = x*0;
% plotting the x and y zero lines (axis)
 plot(x,y,'k');
 plot(y,x,'k');
% axis tight 
% axis equal 
axis([-2 7 -4.5 4.5])

title(" Region of Absolute Stability for Three Step BDF  ");
xlabel("Real Axis");
ylabel("Imaginary");;


hold off