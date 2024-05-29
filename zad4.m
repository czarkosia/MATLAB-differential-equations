clear
close all

H = linspace(1e-4, 1e-2);
P = [13 0.14 0.06 16];
Lotka_Volterra = @(t,u) [u(1).*(P(1)-P(2)*u(2)); u(2).*(P(3)*u(1)-P(4))];
fx = @(x,y) x*(P(1)-P(2)*y);
fy = @(x,y) y*(P(3)*x-P(4)); 

tspan = [0 1];
u0 = [310; 50];
options = odeset('RelTol',1e-8, 'AbsTol', 1e-12);
[tref,u] = ode45(Lotka_Volterra, tspan, u0, options);

error_xa = zeros(length(H), 1);
error_xb = zeros(length(H), 1);
error_xc = zeros(length(H), 1);
error_xd = zeros(length(H), 1);

idx = 1;
for h = H
    clear xa xb xc xd ya yb yc yd xref yref 
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

    %zamknięta metoda Eulera
    for n = 1:N-1
        F = @(u) [xb(n) - u(1) + h*fx(u(1),u(2));
        yb(n) - u(2) + h*fy(u(1),u(2))];
    
        u_temp = fsolve(F, [xb(n); yb(n)]);

        xb(n+1) = u_temp(1);
        yb(n+1) = u_temp(2);
    end

    delta_t = abs(tref-t);
    if (length(tref) > N)
        xref = xa;
        yref = ya;
        for n = 1:N
            [tref,i] = min(delta_t(:,n));
            xref(n) = u(i,1);
            yref(n) = u(i,2);
        end
    else
        xref = u(1);
        yref = u(2);
        ua = [xa ya];
        ub = [xb yb];
        for n = 1:length(tref)
            [t,i] = min(delta_t(n,:));
            xa(n) = ua(i,1);
            ya(n) = ua(i,2);
            xb(n) = ub(i,1);
            yb(n) = ub(i,2);
        end
    end

    error_xa(idx) = sqrt(sum((xa - xref).^2))/sqrt(sum(xref.^2));
    error_xb(idx) = sqrt(sum((xb - xref).^2))/sqrt(sum(xref.^2));
    idx = idx + 1;
end

loglog(H, error_xa, H, error_xb)
