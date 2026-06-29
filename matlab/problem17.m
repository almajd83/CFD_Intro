clear; clc; close all

nx = 201;
L = 1.0;
dx = L/(nx-1);
x = linspace(0,L,nx);

gamma = 1.4;
nt = 200;
CFL = 0.5;

q = zeros(3,nx);

% Initial conditions
for i = 1:nx
    if x(i) < 0.5
        rho = 1.0;
        u = 0.0;
        p = 1.0;
    else
        rho = 0.125;
        u = 0.0;
        p = 0.1;
    end

    E = p/(gamma-1) + 0.5*rho*u^2;

    q(1,i) = rho;
    q(2,i) = rho*u;
    q(3,i) = E;
end

% Time stepping
for n = 1:nt

    rho = q(1,:);
    u = q(2,:)./rho;
    E = q(3,:);
    p = (gamma-1)*(E - 0.5*rho.*u.^2);

    a = sqrt(gamma*p./rho);
    max_speed = max(abs(u) + a);
    dt = CFL * dx / max_speed;

    % Fluxes
    F = zeros(3,nx);
    F(1,:) = rho .* u;
    F(2,:) = rho .* u.^2 + p;
    F(3,:) = u .* (E + p);

    F_half = zeros(3,nx-1);

    % Rusanov (local Lax-Friedrichs) flux
    for i = 1:nx-1
        qL = q(:,i);
        qR = q(:,i+1);

        FL = F(:,i);
        FR = F(:,i+1);

        rhoL = qL(1);
        uL = qL(2)/rhoL;
        pL = (gamma-1)*(qL(3) - 0.5*rhoL*uL^2);
        aL = sqrt(gamma*pL/rhoL);

        rhoR = qR(1);
        uR = qR(2)/rhoR;
        pR = (gamma-1)*(qR(3) - 0.5*rhoR*uR^2);
        aR = sqrt(gamma*pR/rhoR);

        smax = max(abs(uL)+aL, abs(uR)+aR);

        F_half(:,i) = 0.5*(FL + FR) - 0.5*smax*(qR - qL);
    end

    % Update solution
    qn = q;

    for i = 2:nx-1
        q(:,i) = qn(:,i) - dt/dx * (F_half(:,i) - F_half(:,i-1));
    end

    % Boundary conditions (zero gradient)
    q(:,1) = q(:,2);
    q(:,end) = q(:,end-1);
end

% Recover variables
rho = q(1,:);
u = q(2,:)./rho;
E = q(3,:);
p = (gamma-1)*(E - 0.5*rho.*u.^2);

% Plot results
figure

subplot(3,1,1)
plot(x, rho, 'LineWidth', 2)
ylabel('\rho')
title('Density')

subplot(3,1,2)
plot(x, u, 'LineWidth', 2)
ylabel('u')
title('Velocity')

subplot(3,1,3)
plot(x, p, 'LineWidth', 2)
ylabel('p')
xlabel('x')
title('Pressure')