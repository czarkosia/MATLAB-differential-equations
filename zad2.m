clear
close all

h = 0.005;
t = 0:h:1;
P = [13 0.14 0.06 16];

x = zeros(1/h+1, 1);
y = zeros(1/h+1, 1);
x(1) = 310;
y(1) = 50;

%otwarta metoda Eulera
for n = 2:1/h+1
    x(n) = x(n-1) - h*x(n-1).*(P(1)-P(2)*y(n-1));
    y(n) = y(n-1) - h*y(n-1).*(P(3)*x(n-1)-P(4));
end
plot (t, x, t, y)

%zamknięta metoda Eulera