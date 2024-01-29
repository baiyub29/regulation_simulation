% Parameters
L = 1;          % Length of the domain
T = 10;          % Final time
dx = 0.01;
CLF = 0.8;
maxv = 2.5;
dt = CLF * dx / maxv;
Nx = L / dx + 1;       % Number of spatial grid points
Nt = T / dt + 1;       % Number of time steps

% Discretization
x = linspace(0, L, Nx);
t = linspace(0, T, Nt);

% control
s2t = s2(t);
ut = u(t, s2t);

% control deriative
s2_0 = s2t(1);
ud_0 = 1+2.5*(0.5*sin(s2_0)*(cos(0.8*s2_0+0.8)-3*sin(0.5*s2_0+0.5) ...
    -2*sin(0.5*s2_0+0.5)+sin(0.5*s2_0))/(1+0.5*cos(s2_0))^2 ...
    -cos(0.5*s2_0)+(-1.5*cos(0.5*s2_0+0.5)+0.5*cos(0.5*s2_0) ...
    -0.8*sin(0.8*s2_0+0.8)-cos(0.5*s2_0+0.5))/(1+0.5*cos(s2_0)))/(1+cos(2*pi*s2_0));

% Initial conditions for the 2 by 2 system
w = zeros(2, Nx, Nt);

% Set initial condition based on the given boundary condition
w(1, :, 1) = 0.5 * (1 - 3 * ud_0) * x + 1.5 * ut(1) - 0.6 * ud_0 - 2.4;
w(2, :, 1) = 0.4 * (ud_0 - 1) * (x - 1) + ut(1) - 2;

% Discretize the system using finite differences
for n = 1:Nt-1
    % Calculate spatial derivatives
    dw1_dx = diff(w(1, :, n)) / dx;
    dw2_dx = diff(w(2, :, n)) / dx;

    % Update the solution using the provided hyperbolic system
    w(1, 2:Nx, n+1) = w(1, 2:Nx, n) - dt * (1 + 0.5 * sin(2*pi*t(n))) * dw1_dx + dt * cos(0.5*t(n));
    w(2, 1:Nx-1, n+1) = w(2, 1:Nx-1, n) + dt * (1.5 + cos(2*pi*t(n))) * dw2_dx + dt * cos(0.5*t(n));

    % Implement boundary conditions
    w(1, 1, n+1) = (1 + 0.5 * cos(t(n+1))) * w(2, 1, n+1) + sin(0.5*t(n+1));
    w(2, Nx, n+1) = ut(n+1) - 2 * cos(0.5*t(n+1));
end

y = squeeze(w(1, Nx, :))' + 3 * sin(0.5 * t);

% Plot the results
[XX, TT] = meshgrid(x, t);
figure;
s1 = mesh(XX, TT, squeeze(w(1, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$w_1(t,x)$','Interpreter','latex');
%title('Component 1 of the Solution');
%s1.EdgeColor = 'none';
%colorbar;

figure;
s2 = mesh(XX, TT, squeeze(w(2, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$w_2(t,x)$','Interpreter','latex');
%title('Component 2 of the Solution');
%s2.EdgeColor = 'none';

% Additional customization for visualization
%colormap('jet');
%shading interp;
%colorbar;

% plot tracking behaviour
r = cos(0.8 * t);
figure;
subplot(2, 1, 1);
plot(t,y,t,r,'--');
legend('$y(t)$','$r(t)$','Interpreter','latex','FontName','Times New Roman','FontSize',8);

ylabel({'Output $y(t)$';'Reference $r(t)$'},'Interpreter','latex','FontName','Times New Roman','FontSize',8);

% plot tracking error
subplot(2, 1, 2);
plot(t,y-r);

xlabel('Time $t$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);
ylabel('Error $e_y(t)$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);

% plot control function
figure;
plot(t,ut);

xlabel('Time $t$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);
ylabel('Control $u(t)$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);