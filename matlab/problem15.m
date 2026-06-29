div = zeros(ny,nx);

for j = 2:ny-1
    for i = 2:nx-1
        dudx = (u_vel(j,i+1) - u_vel(j,i-1))/(2*dx);
        dvdy = (v_vel(j+1,i) - v_vel(j-1,i))/(2*dy);
        div(j,i) = dudx + dvdy;
    end
end

div_norm = sqrt(sum(sum(div.^2)));

figure
contourf(X,Y,div,30)
colorbar
xlabel('x')
ylabel('y')
title('Divergence Error')

disp(['Divergence norm = ', num2str(div_norm)])