clear
close all

h = 0.005;
N = 1/h + 1;
P = [13 0.14 0.06 16];
fx = @(x,y) x*(P(1)-P(2)*y);
fy = @(x,y) y*(P(3)*x-P(4)); 

t = 0:h:1;
xa = zeros(N, 1);
ya = zeros(N, 1);
t(1) = 0;
xa(1) = 310;
ya(1) = 50;
xb = xa;
yb = ya;
xc = xa;
yc = ya;
xd = xa;
yd = ya;

%otwarta metoda Eulera
for n = 2:N
    xa(n) = xa(n-1) - h*(P(1)*xa(n-1) - P(2)*xa(n-1)*ya(n-1));
    ya(n) = ya(n-1) - h*(P(3)*xa(n-1)*ya(n-1) - P(4)*ya(n-1));
end
figure;
plot (t, xa, t, ya)
legend("x", "y")

%zamknięta metoda Eulera
%powinno się powoli osłabiać czy coś
for n = 2:N
    xb(n) = xb(n-1)/(1 - h*(P(1)-P(2)* yb(n-1)/(1 - h*(P(3)*xb(n)-P(4))) ));
    yb(n) = yb(n-1)/(1 - h*(P(3)*xb(n)-P(4)));
end
figure;
plot(t, xb, t, yb);
legend("x", "y")

%otwarta metoda punktu środkowego
for n = 2:N
    x_temp = xc(n-1) + 1/2*h*xc(n-1)*(P(1) - P(2)*yc(n-1));
    y_temp = yc(n-1) + 1/2*h*yc(n-1)*(P(3)*xc(n-1) - P(4));
    xc(n) = xb(n-1) + h*x_temp;
    yc(n) = yb(n-1) + h*y_temp;
end
figure;
plot(t, xc, t, yc);
legend("x", "y")

%metoda Adamsa-Moultona 3. rzędu
%dwukrokowa, pierwszy xd(2) i yd(2) wyznacz zamkniętą metodą Eulera
xd(2) = (xd(1) + 9/24*h*fx(xd(1), yd(1)) )/( 1 - 19/24*h*(P(3)*xb(n)-P(4)) );
