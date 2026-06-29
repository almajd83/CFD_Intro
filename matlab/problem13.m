clear; clc; close all

nx = 41;
ny = 41;
Lx = 2.0;
Ly = 2.0;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);

nt = 500;
nit = 50;
rho = 1.0;
nu = 0.1;
dt = 0.001;

u_vel = zeros(ny,nx);
v_vel = zeros(ny,nx);
p = zeros(ny,nx);
b = zeros(ny,nx);

for n = 1:nt
    un = u_vel;
    vn = v_vel;

    % Build source term b
    for j = 2:ny-1
        for i = 2:nx-1
            dudx = (un(j,i+1) - un(j,i-1))/(2*dx);
            dudy = (un(j+1,i) - un(j-1,i))/(2*dy);
            dvdx = (vn(j,i+1) - vn(j,i-1))/(2*dx);
            dvdy = (vn(j+1,i) - vn(j-1,i))/(2*dy);

            b(j,i) = rho*((1/dt)*(dudx + dvdy) ...
                      - dudx^2 - 2*dudy*dvdx - dvdy^2);
        end
    end

    % Pressure Poisson
    for k = 1:nit
        pn = p;

        for j = 2:ny-1
            for i = 2:nx-1
                p(j,i) = ((pn(j,i+1) + pn(j,i-1))*dy^2 ...
                        + (pn(j+1,i) + pn(j-1,i))*dx^2 ...
                        - b(j,i)*dx^2*dy^2) ...
                        /(2*(dx^2 + dy^2));
            end
        end

        p(:,1)   = p(:,2);
        p(:,end) = p(:,end-1);
        p(1,:)   = p(2,:);
        p(end,:) = 0.0;
    end

    % Velocity update
    for j = 2:ny-1
        for i = 2:nx-1
            u_vel(j,i) = un(j,i) ...
                - un(j,i)*dt/dx*(un(j,i) - un(j,i-1)) ...
                - vn(j,i)*dt/dy*(un(j,i) - un(j-1,i)) ...
                - dt/(2*rho*dx)*(p(j,i+1) - p(j,i-1)) ...
                + nu*dt/dx^2*(un(j,i+1) - 2*un(j,i) + un(j,i-1)) ...
                + nu*dt/dy^2*(un(j+1,i) - 2*un(j,i) + un(j-1,i));

            v_vel(j,i) = vn(j,i) ...
                - un(j,i)*dt/dx*(vn(j,i) - vn(j,i-1)) ...
                - vn(j,i)*dt/dy*(vn(j,i) - vn(j-1,i)) ...
                - dt/(2*rho*dy)*(p(j+1,i) - p(j-1,i)) ...
                + nu*dt/dx^2*(vn(j,i+1) - 2*vn(j,i) + vn(j,i-1)) ...
                + nu*dt/dy^2*(vn(j+1,i) - 2*vn(j,i) + vn(j-1,i));
        end
    end

    % Boundary conditions
    u_vel(1,:) = 0.0;
    u_vel(:,1) = 0.0;
    u_vel(:,end) = 0.0;
    u_vel(end,:) = 1.0;

    v_vel(1,:) = 0.0;
    v_vel(end,:) = 0.0;
    v_vel(:,1) = 0.0;
    v_vel(:,end) = 0.0;
end

% Plot pressure and velocity field
figure
contourf(X,Y,p,30)
colorbar
hold on
quiver(X,Y,u_vel,v_vel,'k')
xlabel('x')
ylabel('y')
title('Lid-Driven Cavity: Pressure and Velocity')

figure
streamslice(X,Y,u_vel,v_vel)
xlabel('x')
ylabel('y')
title('Lid-Driven Cavity: Streamlines')
