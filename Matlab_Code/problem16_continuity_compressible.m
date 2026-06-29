clear; clc; close all

nx = 101;
L = 1.0;

dx = L/(nx-1);
x = linspace(0,L,nx);

nt = 100;
dt = 0.001;

rho = ones(1,nx);
u = ones(1,nx);

for i = 1:nx
    if x(i) >= 0.4 && x(i) <= 0.6
        rho(i) = 2.0;
    end
end

rho_initial = rho;

for n = 1:nt
    rhon = rho;
    flux = rhon .* u;

    for i = 2:nx
        rho(i) = rhon(i) - dt/dx * (flux(i) - flux(i-1));
    end

    rho(1) = 1.0;
end

plot(x, rho_initial, '--', 'LineWidth', 2)
hold on
plot(x, rho, 'LineWidth', 2)
xlabel('x')
ylabel('\rho')
legend('Initial', 'Final')
title('Compressible Continuity Equation')
grid on