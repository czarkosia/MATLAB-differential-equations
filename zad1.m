clear
close all

P = [13 0.14 0.06 16];
Lotka_Volterra = @(t,u) [u(1).*(P(1)-P(2)*u(2)); u(2).*(P(3)*u(1)-P(4))];

tspan = [0 1];
u0 = [310; 50];
[t,u] = ode45(Lotka_Volterra, tspan, u0);

figure;
plot(t, u);