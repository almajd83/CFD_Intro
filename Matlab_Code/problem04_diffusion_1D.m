clear; clc; close all

nx = 101;
L = 2.0;
dx = L/(nx-1);
x = linspace(0,L,nx);

nt = 50;
nu = 0.05;
dt = 0.25*dx^2/nu;

u_field = ones(1,nx);

for i = 1:nx
    if x(i) >= 0.5 && x(i) <= 1.0
        u_field(i) = 2.0;
    end
end

u_initial = u_field;

for n = 1:nt
    un = u_field;

    for i = 2:nx-1
        u_field(i) = un(i) + nu*dt/dx^2*(un(i+1) - 2*un(i) + un(i-1));
    end

    u_field(1) = 1.0;
    u_field(nx) = 1.0;
end

plot(x,u_initial,'--','LineWidth',2)
hold on
plot(x,u_field,'LineWidth',2)
xlabel('x')
ylabel('u')
legend('Initial','Final')
title('1D Diffusion')
grid on
