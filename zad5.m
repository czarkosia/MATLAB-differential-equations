clear
close all
T = readtable ("MNUB_24L_P3_dane08.csv");

x = T{:,2};
y = T{:,3};
P0 = [13 0.14 0.06 16];
u0 = [x(1), y(1)];

p_est = fminsearch(@(P) Jfun(P, T{:,1}, x, y, u0), P0);

function J = Jfun(P, t, x, y, u0)
    Lotka_Volterra = @(t,u,P) [u(1).*(P(1)-P(2)*u(2)); u(2).*(P(3)*u(1)-P(4))];

    [~, uest] = ode45(@(t, u) Lotka_Volterra(t, u, P), t, u0);
    xest = uest(1);
    yest = uest(2);

    J = sum((xest - x).^2) + sum((yest - y).^2);
end