clear; clc; close all

nx = 81;
ny = 81;
Lx = 2.0;
Ly = 2.0;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);

nt = 100;
nu = 0.05;
dt = 0.0005;

u_field = ones(ny,nx);

for j = 1:ny
    for i = 1:nx
        if x(i) >= 0.5 && x(i) <= 1.0 && y(j) >= 0.5 && y(j) <= 1.0
            u_field(j,i) = 2.0;
        end
    end
end

for n = 1:nt
    un = u_field;

    for j = 2:ny-1
        for i = 2:nx-1
            diffusion_x = (un(j,i+1) - 2*un(j,i) + un(j,i-1)) / dx^2;
            diffusion_y = (un(j+1,i) - 2*un(j,i) + un(j-1,i)) / dy^2;
            u_field(j,i) = un(j,i) + nu*dt*(diffusion_x + diffusion_y);
        end
    end

    u_field(1,:) = 1.0;
    u_field(end,:) = 1.0;
    u_field(:,1) = 1.0;
    u_field(:,end) = 1.0;
end

surf(X,Y,u_field)
xlabel('x')
ylabel('y')
zlabel('u')
title('2D Diffusion')
shading interp