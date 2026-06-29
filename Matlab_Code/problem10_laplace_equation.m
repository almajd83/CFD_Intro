clear; clc; close all

nx = 81;
ny = 81;
Lx = 2.0;
Ly = 1.0;

dx = Lx/(nx-1);
dy = Ly/(ny-1);

x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);

p = zeros(ny,nx);

p(:,1) = 0.0;
p(:,end) = y';
p(1,:) = p(2,:);
p(end,:) = p(end-1,:);

nit = 500;

for k = 1:nit
    pn = p;

    for j = 2:ny-1
        for i = 2:nx-1
            p(j,i) = ((pn(j,i+1) + pn(j,i-1))*dy^2 ...
                    + (pn(j+1,i) + pn(j-1,i))*dx^2) ...
                    / (2*(dx^2 + dy^2));
        end
    end

    p(:,1) = 0.0;
    p(:,end) = y';
    p(1,:) = p(2,:);
    p(end,:) = p(end-1,:);
end

contourf(X,Y,p,20)
colorbar
xlabel('x')
ylabel('y')
title('Laplace Equation Solution')
