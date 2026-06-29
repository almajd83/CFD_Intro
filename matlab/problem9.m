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
dt = 0.001;
nu = 0.01;

u = ones(ny,nx);
v = ones(ny,nx);

for j = 1:ny
    for i = 1:nx
        if x(i) >= 0.5 && x(i) <= 1.0 && y(j) >= 0.5 && y(j) <= 1.0
            u(j,i) = 2.0;
            v(j,i) = 2.0;
        end
    end
end

for n = 1:nt
    un = u;
    vn = v;

    for j = 2:ny-1
        for i = 2:nx-1
            conv_u = un(j,i)*(un(j,i)-un(j,i-1))/dx ...
                   + vn(j,i)*(un(j,i)-un(j-1,i))/dy;

            diff_u = (un(j,i+1)-2*un(j,i)+un(j,i-1))/dx^2 ...
                   + (un(j+1,i)-2*un(j,i)+un(j-1,i))/dy^2;

            conv_v = un(j,i)*(vn(j,i)-vn(j,i-1))/dx ...
                   + vn(j,i)*(vn(j,i)-vn(j-1,i))/dy;

            diff_v = (vn(j,i+1)-2*vn(j,i)+vn(j,i-1))/dx^2 ...
                   + (vn(j+1,i)-2*vn(j,i)+vn(j-1,i))/dy^2;

            u(j,i) = un(j,i) - dt*conv_u + nu*dt*diff_u;
            v(j,i) = vn(j,i) - dt*conv_v + nu*dt*diff_v;
        end
    end

    u(1,:) = 1.0; u(end,:) = 1.0; u(:,1) = 1.0; u(:,end) = 1.0;
    v(1,:) = 1.0; v(end,:) = 1.0; v(:,1) = 1.0; v(:,end) = 1.0;
end

figure
contourf(X,Y,u,20)
colorbar
title('2D Burgers: u component')
xlabel('x')
ylabel('y')

figure
contourf(X,Y,v,20)
colorbar
title('2D Burgers: v component')
xlabel('x')
ylabel('y')
