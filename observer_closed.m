% Parameters
T = 10; % simulation time
Nx = 100; % Number of spatial grid points
Nt = 3000; % Number of time steps

% Discretization
x = linspace(0, 1, Nx); 
t = linspace(0, T, Nt); 
dx = x(2) - x(1);
dt = t(2) - t(1);

A = [0, -0.5; 0.5, 0];

% CFL condition
lambda1 = max(abs(1 + 0.5 * sin(2 * pi * t))); % speed 1
lambda2 = max(abs(1.5 + cos(2 * pi * t))); % speed 2
CFL = max(lambda1, lambda2) * dt / dx;
if CFL >= 1
    error('CFL condition violated, reduce dt or increase Nx.');
end

[tt, xx] = meshgrid(t, x);

s1intx = zeros(Nx, Nt);
s2intx = zeros(Nx, Nt);
for n = 1:Nx
    s1intx(n,:) = s1in_tx(t,x(n)*ones(1, Nt));
    s2intx(n,:) = s2in_tx(t,x(n)*ones(1, Nt));
end

% calculate \hat{N}(t,x)
hatN_11 = 2*sin(0.5*tt)-sin(0.5*s1intx);
hatN_12 = 2*cos(0.5*tt)-cos(0.5*s1intx);
hatN_21 = 2*sin(0.5*tt)-2*cos(0.5*s2intx)-sin(0.5*s2intx-0.5);
hatN_22 = 2*cos(0.5*tt)+2*sin(0.5*s2intx)-cos(0.5*s2intx-0.5);

s2int0 = s2intx(1,:);

% check the observability condition
dst = (1.5+cos(2*pi*t))./(1.5+cos(2*pi*s2int0));
dhatN_21 = cos(0.5*t)+dst.*(sin(0.5*s2int0)-0.5*cos(0.5*s2int0-0.5));
dhatN_22 = -sin(0.5*t)+dst.*(cos(0.5*s2int0)+0.5*sin(0.5*s2int0-0.5));

detN = hatN_21(1,:).*dhatN_22 - hatN_22(1,:).*dhatN_21;

if all(abs(detN) > 0.05)
    disp('The system is observable');
else
    error('The system is not observable. Terminating the program.');
end

% calculate observer gains ld1, ld2, and lw(t,x)=hat{N}(t,x)*exp(-t*S_d)*ld1
N_1 = hatN_21(1,:).*cos(0.5*t) - hatN_22(1,:).*sin(0.5*t);
N_2 = hatN_21(1,:).*sin(0.5*t) + hatN_22(1,:).*cos(0.5*t);

%tN_1 = -2*cos(0.5*(t-s2int0))+sin(0.5*(t-s2int0+1));
%tN_2 = 2-2*sin(0.5*(t-s2int0))-cos(0.5*(t-s2int0+1));

ld1 = zeros(2, Nt);
ld2 = zeros(2, Nt);

for n=1:Nt
    ld1(:,n) = place(A', [N_1(n), N_2(n)]', [-1+1i, -1-1i])';
    ld2(:,n) = place(A', [N_1(n), N_2(n)]', [-2+1i, -2-1i])';
end

%%
% true disturbance system vd=(cos(0.5t),sin(0.5t)),vd0=(1,0)
% delay D_d = 1; t(n) ---- t(n-300)

% s_2^out(t,1)
s2t = s2(t);

% control
ut = zeros(1,Nt);

% state for the 2 by 2 system
w = zeros(2, Nx, Nt);

% state for the observer system
hat_vd = zeros(2, Nt);
hat_mud = zeros(2, Nt);
hat_w = zeros(2, Nx, Nt);

% state for the finite time observer system
hat_vdf = zeros(2, Nt);
hat_wf = zeros(2, Nx, Nt);

% Initial conditions for the observer ode system
hat_vd(:, 1) = [0; 0];
%hat_mud(:, 1) = [2; -1];
hat_mud(:, 1) = [0; 0];

Phi_10 = expm((A-ld1(:,1)*[N_1(1),N_2(1)])*300*dt);
Phi_20 = expm((A-ld2(:,1)*[N_1(1),N_2(1)])*300*dt);

hat_vdf(:, 1) = [eye(2), zeros(2)] * ...
            ([eye(2), Phi_10; eye(2), Phi_20] \ ...
            [hat_vd(:,1)-Phi_10*hat_vd(:,1); hat_mud(:,1)-Phi_20*hat_mud(:,1)]);

ini_v_in_u = 2*(1-sin(0.5*s2t(1)))*hat_vdf(1,1)+2*(1-cos(0.5*s2t(1)))*hat_vdf(2,1) +...
    (cos(0.8*(s2t(1)+1)) + (sin(0.5*s2t(1))-5*sin(0.5*(s2t(1)+1)))*hat_vdf(1,1) + ...
    (cos(0.5*s2t(1))-5*cos(0.5*(s2t(1)+1)))*hat_vdf(2,1)) / (1+0.5*cos(s2t(1)));

wf_in_w_1 = hatN_11(:,1)*(hat_vdf(1, 1) - hat_vd(1, 1))+...
    hatN_12(:,1)*(hat_vdf(2, 1) - hat_vd(2, 1));
wf_in_w_2 = hatN_21(:,1)*(hat_vdf(1, 1) - hat_vd(1, 1))+...
    hatN_22(:,1)*(hat_vdf(2, 1) - hat_vd(2, 1));

% Initial conditions for the 2 by 2 system are 0
% Initial conditions for the observer pde system 

hat_w(1, :, 1) = (-hat_vd(2, 1)-wf_in_w_1(Nx)+ini_v_in_u-2) * x + hat_vd(2, 1);
hat_w(2, :, 1) = (-wf_in_w_1(Nx)+ini_v_in_u- 2* hat_vd(1, 1)) * x;

hat_wf(1, :, 1) = hat_w(1, :, 1) + wf_in_w_1';
hat_wf(2, :, 1) = hat_w(2, :, 1) + wf_in_w_2';

ut(1) = -hat_wf(1, Nx, 1)+ini_v_in_u;



% Discretize the system using finite differences
for n = 1:Nt-1
    Phi_1 = eye(2);
    Phi_2 = eye(2);
    for k = (n-299):(n)
        if k >= 1
            A_k = A ...
                -[ld1(1,k);ld1(2,k)] * [N_1(k), N_2(k)];
            B_k = [0, -0.5; 0.5, 0] ...
                -[ld2(1,k);ld2(2,k)] * [N_1(k), N_2(k)];
        else
            A_k = A ...
                -[ld1(1,1);ld1(2,1)] * [N_1(1), N_2(1)];
            B_k = [0, -0.5; 0.5, 0] ...
                -[ld2(1,1);ld2(2,1)] * [N_1(1), N_2(1)];
        end
        Phi_1k = expm(A_k * dt);  
        Phi_1 = Phi_1k * Phi_1;  
        Phi_2k = expm(B_k * dt);  
        Phi_2 = Phi_2k * Phi_2; 
    end

    % Calculate spatial derivatives
    dw1_dx = diff(w(1, :, n)) / dx;
    dw2_dx = diff(w(2, :, n)) / dx;

    dhatw1_dx = diff(hat_w(1, :, n)) / dx;
    dhatw2_dx = diff(hat_w(2, :, n)) / dx;

    % Update the solution using the provided hyperbolic system
    w(1, 2:Nx, n+1) = w(1, 2:Nx, n) - dt * (1 + 0.5 * sin(2*pi*t(n))) * dw1_dx + dt * cos(0.5*t(n));
    w(2, 1:Nx-1, n+1) = w(2, 1:Nx-1, n) + dt * (1.5 + cos(2*pi*t(n))) * dw2_dx + dt * cos(0.5*t(n));

    hat_w(1, 2:Nx, n+1) = hat_w(1, 2:Nx, n) - dt * (1 + 0.5 * sin(2*pi*t(n))) * dhatw1_dx ...
        + dt * (hat_vd(1, n) * ones(1,Nx-1)...
        + ((cos(0.5*t(n))*ld1(1,n) + sin(0.5*t(n))*ld1(2,n))*hatN_11(2:Nx,n).'...
        + (-sin(0.5*t(n))*ld1(1,n) + cos(0.5*t(n))*ld1(2,n))*hatN_12(2:Nx,n).') ...
        * (w(2, 1, n) - hat_w(2, 1, n)));
    hat_w(2, 1:Nx-1, n+1) = hat_w(2, 1:Nx-1, n) + dt * (1.5 + cos(2*pi*t(n))) * dhatw2_dx ...
        + dt * (hat_vd(1, n) * ones(1,Nx-1)...
        + ((cos(0.5*t(n))*ld1(1,n) + sin(0.5*t(n))*ld1(2,n))*hatN_21(1:Nx-1,n).'...
        + (-sin(0.5*t(n))*ld1(1,n) + cos(0.5*t(n))*ld1(2,n))*hatN_22(1:Nx-1,n).') ...
        * (w(2, 1, n) - hat_w(2, 1, n)));

    hat_vd(1, n+1) = hat_vd(1, n) + dt * (-0.5 * hat_vd(2, n) ...
        + ld1(1, n) * (w(2, 1, n) - hat_w(2, 1, n)));
    hat_vd(2, n+1) = hat_vd(2, n) + dt * (0.5 * hat_vd(1, n) ...
        + ld1(2, n) * (w(2, 1, n) - hat_w(2, 1, n)));
    hat_mud(1, n+1) = hat_mud(1, n) + dt * (-0.5 * hat_mud(2, n) ...
        + ld2(1, n) * (w(2, 1, n) - hat_w(2, 1, n)...
        + N_1(n) * (hat_vd(1,n) - hat_mud(1, n))...
        + N_2(n) * (hat_vd(2,n) - hat_mud(2, n))));
    hat_mud(2, n+1) = hat_mud(2, n) + dt * (0.5 * hat_mud(1, n) ...
        + ld2(2, n) * (w(2, 1, n) - hat_w(2, 1, n)...
        + N_1(n) * (hat_vd(1,n) - hat_mud(1, n))...
        + N_2(n) * (hat_vd(2,n) - hat_mud(2, n))));

    % calculate the finite time ode observer and control
    if n >= 300
        hat_vdf(:, n+1) = [eye(2), zeros(2)] * ...
            ([eye(2), Phi_1; eye(2), Phi_2] \ ...
            [hat_vd(:,n+1)-Phi_1*hat_vd(:,n-299); hat_mud(:,n+1)-Phi_2*hat_mud(:,n-299)]);
    else
        hat_vdf(:, n+1) = [eye(2), zeros(2)] * ...
            ([eye(2), Phi_1; eye(2), Phi_2] \ ...
            [hat_vd(:,n+1)-Phi_1*hat_vd(:,1); hat_mud(:,n+1)-Phi_2*hat_mud(:,1)]);
    end

    esthatw1_t1 = hat_w(1, Nx, n+1) + ...
        (hatN_11(Nx,n+1)*cos(0.5*t(n+1))-hatN_12(Nx,n+1)*sin(0.5*t(n+1)))*...
        (hat_vdf(1, n+1) - hat_vd(1, n+1)) + ...
        (hatN_11(Nx,n+1)*sin(0.5*t(n+1))+hatN_12(Nx,n+1)*cos(0.5*t(n+1)))*...
        (hat_vdf(2, n+1) - hat_vd(2, n+1));

    ut(n+1) = u_s(t(n+1), s2t(n+1), esthatw1_t1, hat_vdf(1, n+1), hat_vdf(2, n+1));

    % Implement boundary conditions
    w(1, 1, n+1) = (1 + 0.5 * cos(t(n+1))) * w(2, 1, n+1) + sin(0.5*t(n+1));
    w(2, Nx, n+1) = w(1, Nx, n+1) + ut(n+1) - 2 * cos(0.5*t(n+1));

    hat_w(1, 1, n+1) = (1 + 0.5 * cos(t(n+1))) * w(2, 1, n+1) + hat_vd(2, n+1);
    hat_w(2, Nx, n+1) = hat_w(1, Nx, n+1) + ut(n+1) - 2 * hat_vd(1, n+1);

    % finite time pde observer
    hat_wf(1, :, n+1) = hat_w(1, :, n+1) + ...
        (hatN_11(:,n+1)'*cos(0.5*t(n+1))-hatN_12(:,n+1)'*sin(0.5*t(n+1)))*...
        (hat_vdf(1, n+1) - hat_vd(1, n+1)) + ...
        (hatN_11(:,n+1)'*sin(0.5*t(n+1))+hatN_12(:,n+1)'*cos(0.5*t(n+1)))*...
        (hat_vdf(2, n+1) - hat_vd(2, n+1));
    hat_wf(2, :, n+1) = hat_w(2, :, n+1) + ...
        (hatN_21(:,n+1)'*cos(0.5*t(n+1))-hatN_22(:,n+1)'*sin(0.5*t(n+1)))*...
        (hat_vdf(1, n+1) - hat_vd(1, n+1)) + ...
        (hatN_21(:,n+1)'*sin(0.5*t(n+1))+hatN_22(:,n+1)'*cos(0.5*t(n+1)))*...
        (hat_vdf(2, n+1) - hat_vd(2, n+1));
end

y = squeeze(w(1, Nx, :))' + 3 * sin(0.5 * t);


% Plot the results
[XX, TT] = meshgrid(x, t);
figure;
f1 = mesh(XX, TT, squeeze(w(1, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$w_1(t,x)$','Interpreter','latex');

figure;
f2 = mesh(XX, TT, squeeze(hat_wf(1, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$\hat{w}^+_1(t,x)$','Interpreter','latex');

figure;
f3 = mesh(XX, TT, squeeze(w(2, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$w_2(t,x)$','Interpreter','latex');

figure;
f4 = mesh(XX, TT, squeeze(hat_wf(2, :, :))');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
zlabel('$\hat{w}^+_2(t,x)$','Interpreter','latex');

% plot disturbance behaviour
vd1 = cos(0.5 * t);
vd2 = sin(0.5 * t);
figure;
subplot(2, 1, 1);
plot(t,hat_vdf(1,:),t,vd1,'--');
legend('$\hat{v}^+_{d,1}(t)$','$v_{d,1}(t)$','Interpreter','latex','FontName','Times New Roman','FontSize',8);
ylabel({'Disturbance $v_{d,1}(t)$';'Observer $\hat{v}^+_{d,1}(t)$'},'Interpreter','latex','FontName','Times New Roman','FontSize',8);

subplot(2, 1, 2);
plot(t,hat_vdf(1,:)-vd1);

xlabel('Time $t$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);
ylabel('Error $\hat{v}^+_{d,1}(t)-v_{d,1}(t)$','Interpreter','latex','FontName','Times New Roman','FontSize',8);

figure;
subplot(2, 1, 1);
plot(t,hat_vdf(2,:),t,vd2,'--');
legend('$\hat{v}^+_{d,2}(t)$','$v_{d,2}(t)$','Interpreter','latex','FontName','Times New Roman','FontSize',8);
ylabel({'Disturbance $v_{d,2}(t)$';'Observer $\hat{v}^+_{d,2}(t)$'},'Interpreter','latex','FontName','Times New Roman','FontSize',8);

subplot(2, 1, 2);
plot(t,hat_vdf(2,:)-vd2);

xlabel('Time $t$','Interpreter','latex', 'FontName','Times New Roman','FontSize',8);
ylabel('Error $\hat{v}^+_{d,2}(t)-v_{d,2}(t)$','Interpreter','latex','FontName','Times New Roman','FontSize',8);

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