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

nt = 50;
dt = 0.001;

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

    for j = 2:ny
        for i = 2:nx
            u(j,i) = un(j,i) ...
                - un(j,i)*dt/dx*(un(j,i)-un(j,i-1)) ...
                - vn(j,i)*dt/dy*(un(j,i)-un(j-1,i));

            v(j,i) = vn(j,i) ...
                - un(j,i)*dt/dx*(vn(j,i)-vn(j,i-1)) ...
                - vn(j,i)*dt/dy*(vn(j,i)-vn(j-1,i));
        end
    end

    u(1,:) = 1.0; u(end,:) = 1.0; u(:,1) = 1.0; u(:,end) = 1.0;
    v(1,:) = 1.0; v(end,:) = 1.0; v(:,1) = 1.0; v(:,end) = 1.0;
end

figure
surf(X,Y,u)
title('u velocity')
xlabel('x')
ylabel('y')
zlabel('u')
shading interp

figure
surf(X,Y,v)
title('v velocity')
xlabel('x')
ylabel('y')
zlabel('v')
shading interp