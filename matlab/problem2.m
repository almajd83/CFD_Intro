clear; clc; close all

nx = 101;
L = 2.0;
dx = L/(nx-1);
x = linspace(0,L,nx);

c = 1.0;
nt = 50;

sigma_values = [0.2, 0.5, 1.0, 1.2];

figure

for case_id = 1:length(sigma_values)

    sigma = sigma_values(case_id);
    dt = sigma*dx/c;

    u = ones(1,nx);

    for i = 1:nx
        if x(i) >= 0.5 && x(i) <= 1.0
            u(i) = 2.0;
        end
    end

    for n = 1:nt
        un = u;

        for i = 2:nx
            u(i) = un(i) - sigma*(un(i) - un(i-1));
        end

        u(1) = 1.0;
    end

    subplot(2,2,case_id)
    plot(x,u,'LineWidth',2)
    title(['CFL = ', num2str(sigma)])
    xlabel('x')
    ylabel('u')
    grid on

end