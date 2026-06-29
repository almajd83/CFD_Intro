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

rho = 1.0;
dt = 0.001;

u = sin(pi*X).*sin(pi*Y);
v = cos(pi*X).*cos(pi*Y);
p = zeros(ny,nx);
b = zeros(ny,nx);

for j = 2:ny-1
    for i = 2:nx-1
        dudx = (u(j,i+1) - u(j,i-1)) / (2*dx);
        dudy = (u(j+1,i) - u(j-1,i)) / (2*dy);
        dvdx = (v(j,i+1) - v(j,i-1)) / (2*dx);
        dvdy = (v(j+1,i) - v(j-1,i)) / (2*dy);

        b(j,i) = rho*((dudx + dvdy)/dt ...
               - dudx^2 - 2*dudy*dvdx - dvdy^2);
    end
end

nit = 500;

for k = 1:nit
    pn = p;

    for j = 2:ny-1
        for i = 2:nx-1
            p(j,i) = ((pn(j,i+1) + pn(j,i-1))*dy^2 ...
                    + (pn(j+1,i) + pn(j-1,i))*dx^2 ...
                    - b(j,i)*dx^2*dy^2) ...
                    / (2*(dx^2 + dy^2));
        end
    end

    p(:,1) = p(:,2);
    p(:,end) = p(:,end-1);
    p(1,:) = p(2,:);
    p(end,:) = 0.0;
end

contourf(X,Y,p,30)
colorbar
xlabel('x')
ylabel('y')
title('Pressure Poisson Equation')