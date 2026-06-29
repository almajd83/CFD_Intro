clear; clc; close all

nx = 101;
L = 2.0;
dx = L/(nx-1);
x = linspace(0,L,nx);

nt = 50;
dt = 0.005;

u = ones(1,nx);

for i = 1:nx
    if x(i) >= 0.5 && x(i) <= 1.0
        u(i) = 2.0;
    end
end

u_initial = u;

for n = 1:nt
    un = u;

    for i = 2:nx
        u(i) = un(i) - un(i)*dt/dx*(un(i) - un(i-1));
    end

    u(1) = 1.0;
end

plot(x,u_initial,'--','LineWidth',2)
hold on
plot(x,u,'LineWidth',2)
xlabel('x')
ylabel('u')
legend('Initial','Final')
title('1D Nonlinear Convection')
grid on
