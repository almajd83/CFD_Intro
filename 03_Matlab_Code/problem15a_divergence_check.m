clear; clc; close all

% Reynolds study
Re_values = [10, 100, 400];
U = 1.0;
L = 2.0;

% Grid
nx = 41;
ny = 41;

dx = L/(nx-1);
dy = L/(ny-1);

x = linspace(0,L,nx);
y = linspace(0,L,ny);
[X,Y] = meshgrid(x,y);

nt = 500;
nit = 50;
rho = 1.0;
dt = 0.001;

divergence_history = zeros(length(Re_values),1);

for case_id = 1:length(Re_values)

    Re = Re_values(case_id);
    nu = U*L/Re;

    % Initialize fields
    u_vel = zeros(ny,nx);
    v_vel = zeros(ny,nx);
    p = zeros(ny,nx);
    b = zeros(ny,nx);

    for n = 1:nt
        un = u_vel;
        vn = v_vel;

        % Build source term
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

        % Pressure Poisson solver
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

            % Pressure BCs
            p(:,1)   = p(:,2);
            p(:,end) = p(:,end-1);
            p(1,:)   = p(2,:);
            p(end,:) = 0.0;
        end

        % Velocity update
        for j = 2:ny-1
            for i = 2:nx-1
                u_vel(j,i) = un(j,i) ...
                    - un(j,i)*dt/dx*(un(j,i)-un(j,i-1)) ...
                    - vn(j,i)*dt/dy*(un(j,i)-un(j-1,i)) ...
                    - dt/(2*rho*dx)*(p(j,i+1)-p(j,i-1)) ...
                    + nu*dt/dx^2*(un(j,i+1)-2*un(j,i)+un(j,i-1)) ...
                    + nu*dt/dy^2*(un(j+1,i)-2*un(j,i)+un(j-1,i));

                v_vel(j,i) = vn(j,i) ...
                    - un(j,i)*dt/dx*(vn(j,i)-vn(j,i-1)) ...
                    - vn(j,i)*dt/dy*(vn(j,i)-vn(j-1,i)) ...
                    - dt/(2*rho*dy)*(p(j+1,i)-p(j-1,i)) ...
                    + nu*dt/dx^2*(vn(j,i+1)-2*vn(j,i)+vn(j,i-1)) ...
                    + nu*dt/dy^2*(vn(j+1,i)-2*vn(j,i)+vn(j-1,i));
            end
        end

        % Boundary conditions (lid-driven cavity)
        u_vel(1,:) = 0.0;
        u_vel(:,1) = 0.0;
        u_vel(:,end) = 0.0;
        u_vel(end,:) = U;

        v_vel(1,:) = 0.0;
        v_vel(end,:) = 0.0;
        v_vel(:,1) = 0.0;
        v_vel(:,end) = 0.0;
    end

    % ---------- Plot results ----------
    figure
    contourf(X,Y,p,30)
    colorbar
    hold on
    quiver(X,Y,u_vel,v_vel,'k')
    title(['Re = ', num2str(Re), ' Pressure & Velocity'])
    xlabel('x')
    ylabel('y')

    figure
    streamslice(X,Y,u_vel,v_vel)
    title(['Re = ', num2str(Re), ' Streamlines'])
    xlabel('x')
    ylabel('y')

    % Velocity magnitude
    vel_mag = sqrt(u_vel.^2 + v_vel.^2);

    figure
    contourf(X,Y,vel_mag,30)
    colorbar
    title(['Re = ', num2str(Re), ' Velocity Magnitude'])
    xlabel('x')
    ylabel('y')

    % ---------- Divergence check ----------
    div = zeros(ny,nx);

    for j = 2:ny-1
        for i = 2:nx-1
            dudx = (u_vel(j,i+1) - u_vel(j,i-1))/(2*dx);
            dvdy = (v_vel(j+1,i) - v_vel(j-1,i))/(2*dy);
            div(j,i) = dudx + dvdy;
        end
    end

    div_norm = sqrt(sum(sum(div.^2)));
    divergence_history(case_id) = div_norm;

    figure
    contourf(X,Y,div,30)
    colorbar
    title(['Re = ', num2str(Re), ' Divergence Error'])
    xlabel('x')
    ylabel('y')

    disp(['Re = ', num2str(Re), ...
          ' | Divergence norm = ', num2str(div_norm)])

end

% ---------- Comparison plot ----------
figure
plot(Re_values, divergence_history,'-o','LineWidth',2)
xlabel('Reynolds number')
ylabel('Divergence norm')
title('Incompressibility error vs Reynolds number')
grid on
