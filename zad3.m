clear
run ("zad2.m")
close all

P = [13 0.14 0.06 16];
Lotka_Volterra = @(t,u) [u(1).*(P(1)-P(2)*u(2)); u(2).*(P(3)*u(1)-P(4))];

u0 = [310; 50];
options = odeset('RelTol',1e-8, 'AbsTol', 1e-12);
[tref,u] = ode45(Lotka_Volterra, t, u0, options);
xref = u(:,1);
yref = u(:,2);

error_xa = sqrt(sum((xa - xref).^2))/sqrt(sum(xref.^2));
error_xb = sqrt(sum((xb - xref).^2))/sqrt(sum(xref.^2));
error_xc = sqrt(sum((xc - xref).^2))/sqrt(sum(xref.^2));
error_xd = sqrt(sum((xd - xref).^2))/sqrt(sum(xref.^2));