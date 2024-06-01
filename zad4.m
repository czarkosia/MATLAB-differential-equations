clear
close all

H = logspace(-4, -2);
P = [13 0.14 0.06 16];

Lotka_Volterra = @(t,u) [u(1).*(P(1)-P(2)*u(2)); u(2).*(P(3)*u(1)-P(4))];
fx = @(x,y) x*(P(1)-P(2)*y);
fy = @(x,y) y*(P(3)*x-P(4)); 
options = odeset('RelTol',1e-8, 'AbsTol', 1e-12);
u0 = [310; 50];

error_xa = zeros(length(H), 1);
error_xb = zeros(length(H), 1);
error_xc = zeros(length(H), 1);
error_xd = zeros(length(H), 1);

idx = 1;
for h = H
    t = 0:h:1;
    N = length(t);
    xa = zeros(N, 1);
    ya = zeros(N, 1);
    xa(1) = u0(1);
    ya(1) = u0(2);
    xb = xa;
    yb = ya;
    xc = xa;
    yc = ya;
    xd = xa;
    yd = ya;

    [tref,u] = ode45(Lotka_Volterra, t, u0, options);
    xref = u(:,1);
    yref = u(:,2);

    %otwarta metoda Eulera
    for n = 2:N
        xa(n) = xa(n-1) + h*fx(xa(n-1),ya(n-1));
        ya(n) = ya(n-1) + h*fy(xa(n-1),ya(n-1));
    end

    %zamknięta metoda Eulera
    for n = 1:N-1
        F = @(uf) [xb(n) - uf(1) + h*fx(uf(1),uf(2));
            yb(n) - uf(2) + h*fy(uf(1),uf(2))];
        
        u_temp = fsolve(F, [xb(n); yb(n)]);
    
        xb(n+1) = u_temp(1);
        yb(n+1) = u_temp(2);
    end

    %otwarta metoda punktu środkowego
    for n = 2:N
        x_temp = xc(n-1) + 1/2*h*fx(xc(n-1),yc(n-1));
        y_temp = yc(n-1) + 1/2*h*fy(xc(n-1),yc(n-1));
        xc(n) = xc(n-1) + h*fx(x_temp,y_temp);
        yc(n) = yc(n-1) + h*fy(x_temp,y_temp);
    end

    %metoda Adamsa-Moultona 3. rzędu
    F = @(uf) [xd(1) - uf(1) + h*fx(uf(1),uf(2));
            yd(1) - uf(2) + h*fy(uf(1),uf(2))];
    u_temp = fsolve(F, [xd(1); yd(1)]);
    xd(2) = u_temp(1);
    yd(2) = u_temp(2);
    
    for n = 2:N-1
        F = @(uf) [xd(n) - uf(1) + 5/12*h*fx(uf(1),uf(2)) + 2/3*h*fx(xd(n),yd(n)) - 1/12*h*fx(xd(n-1),yd(n-1));
            yd(n) - uf(2) + 5/12*h*fy(uf(1),uf(2)) + 2/3*h*fy(xd(n),yd(n)) - 1/12*h*fy(xd(n-1),yd(n-1))];
    
        u_temp = fsolve(F, [xd(n); yd(n)]);
    
        xd(n+1) = u_temp(1);
        yd(n+1) = u_temp(2);
    end

    %zagregowany błąd względny x dla każdej metody
    error_xa(idx) = sqrt(sum((xa - xref).^2))/sqrt(sum(xref.^2));
    error_xb(idx) = sqrt(sum((xb - xref).^2))/sqrt(sum(xref.^2));
    error_xc(idx) = sqrt(sum((xc - xref).^2))/sqrt(sum(xref.^2));
    error_xd(idx) = sqrt(sum((xd - xref).^2))/sqrt(sum(xref.^2));
    idx = idx + 1;
end

loglog(H, error_xa, H, error_xb, H, error_xc, H, error_xd)
legend("Otwarta metoda Eulera", "Zamknięta metoda Eulera", "Otwarta metoda punktu środkowego", "Metoda Adamsa-Moultona 3. rzędu")
xlabel("Wartość kroku h")
ylabel("Zagregowany błąd względny x")