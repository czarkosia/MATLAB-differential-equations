clear
close all

h = 0.005;
P = [13 0.14 0.06 16];
fx = @(x,y) x*(P(1)-P(2)*y);
fy = @(x,y) y*(P(3)*x-P(4)); 

t = 0:h:1;
N = length(t);
xa = zeros(N, 1);
ya = zeros(N, 1);
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
    xa(n) = xa(n-1) + h*fx(xa(n-1),ya(n-1));
    ya(n) = ya(n-1) + h*fy(xa(n-1),ya(n-1));
end
figure;
plot (t, xa, t, ya)
legend("x", "y")

%zamknięta metoda Eulera
for n = 1:N-1
    F = @(u) [xb(n) - u(1) + h*fx(u(1),u(2));
        yb(n) - u(2) + h*fy(u(1),u(2))];
    
    u_temp = fsolve(F, [xb(n); yb(n)]);

    xb(n+1) = u_temp(1);
    yb(n+1) = u_temp(2);
end
figure;
plot(t, xb, t, yb);
legend("x", "y")

%otwarta metoda punktu środkowego
for n = 2:N
    x_temp = xc(n-1) + 1/2*h*fx(xc(n-1),yc(n-1));
    y_temp = yc(n-1) + 1/2*h*fy(xc(n-1),yc(n-1));
    xc(n) = xc(n-1) + h*fx(x_temp,y_temp);
    yc(n) = yc(n-1) + h*fy(x_temp,y_temp);
end
figure;
plot(t, xc, t, yc);
legend("x", "y")

%metoda Adamsa-Moultona 3. rzędu
%dwukrokowa, pierwszy xd(2) i yd(2) wyznacz zamkniętą metodą Eulera
F = @(u) [xd(1) - u(1) + h*fx(u(1),u(2));
        yd(1) - u(2) + h*fy(u(1),u(2))];
u_temp = fsolve(F, [xd(1); yd(1)]);

xd(2) = u_temp(1);
yd(2) = u_temp(2);

% F = @(u) [xd(2) - u(1) + 1/2*h*fx(u(1),u(2)) + 1/2*h*fx(xd(2),yd(2));
%         yd(2) - u(2) + h*fy(u(1),u(2)) + 1/2*h*fy(xd(2),yd(2))];
% u_temp = fsolve(F, [xd(1); yd(1)]);
% 
% xd(3) = u_temp(1);
% yd(3) = u_temp(2);

for n = 2:N-1
    F = @(u) [xd(n) - u(1) + 5/12*h*fx(u(1),u(2)) + 2/3*h*fx(xd(n),yd(n)) - 1/12*h*fx(xd(n-1),yd(n-1));
        yd(n) - u(2) + 5/12*h*fy(u(1),u(2)) + 2/3*h*fy(xd(n),yd(n)) - 1/12*h*fy(xd(n-1),yd(n-1))];

    u_temp = fsolve(F, [xd(n); yd(n)]);

    xd(n+1) = u_temp(1);
    yd(n+1) = u_temp(2);
end
figure;
plot(t, xd, t, yd);
legend("x", "y")

%zad4 najlepiej loglogiem wyrysować
