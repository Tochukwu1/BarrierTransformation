close all; clear all; clc;

% Simulation Parameters
P.Q = diag([1, 1, 1, 1]);
P.R = diag([1, 1]);
P.N = 100;
P.l = 10; 
P.n = 4;
P.q = 4; % number of unknown parameters
P.kc1 = 0.1;
P.kc2 = 10;
P.ka1 = 20;
P.ka2 = 0.2;
P.Beta = 0.8;
P.nu = 100;
P.Beta_1 = diag([100, 100, 100, 100]);
theta = [5.3;1.1;8.45;2.35];

P.p1 = 3.473;
P.p2 = 0.196;
P.p3 = 0.242;

% Initial Conditions
x_0 = [-5;-5;5;5];
s_0 =  [b(x_0(1),-7,5);b(x_0(2),-7,5);b(x_0(3),-5,7);b(x_0(4),-5,7)];
theta_hat0 = [5;5;5;5];
gamma_0 = 10*eye(10);
wa_hat0 = [60;2;2;2;2;2;40;2;2;2];
wc_hat0 = wa_hat0;
Y_x0 = zeros(16,1);
Yf0 = zeros(16,1);
Gs0 = zeros(4,1);
Xf0 = zeros(4,1);

y0 = [s_0;wa_hat0;wc_hat0;gamma_0(:);theta_hat0(:);Y_x0;Yf0(:);Gs0(:);Xf0;x_0];
tspan = [0 5];

% Simulation

options = odeset('OutputFcn',@odeplot, 'OutputSel', P.n+P.l+P.l+P.l^2+1:P.n+P.l+P.l+P.l^2+P.q);
[t,y] = ode45(@(t,y) closeLoopDynamics(t,y, P), tspan, y0, options);

% Plots
plot(t, y(:, P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+P.n+P.q+1:end));
legend('$x_1$', '$x_2$', 'Interpreter','Latex')
title('Graph of x', 'Interpreter','Latex')
grid on

figure
plot(t, y(:, 1:P.n));
legend('$s_1$', '$s_2$', 'Interpreter','Latex')
title('Graph of s', 'Interpreter','Latex')
grid on

figure
plot(t, y(:, P.n+P.l+P.l+P.l^2+1:P.n+P.l+P.l+P.l^2+P.q))
yline(theta(1), '--','Color', 'cyan')
yline(theta(2), '--','Color', 'green')
yline(theta(3), '--','Color', 'magenta')
yline(theta(4), '--','Color', 'blue')
title('Graph of $\theta$', 'Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
ylabel('$\theta$ (rad)','Interpreter','Latex')
legend('$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','Interpreter','Latex');
grid on 

figure
plot(t, y(:, P.n+1:P.n+P.l));
title('Policy Weights', 'FontSize',20,'Interpreter','latex')
xlabel('$Time (s)$','FontSize',20,'Interpreter','latex')
ylabel('$\hat{W}_{a}(t)$','FontSize',20,'Interpreter','latex')
legend('$\hat{W}_{a,1}(t)$','$\hat{W}_{a,2}(t)$', '$\hat{W}_{a,3}(t)$',...
    '$\hat{W}_{a,4}(t)$','$\hat{W}_{a,5}(t)$', '$\hat{W}_{a,6}(t)$',...
    '$\hat{W}_{a,7}(t)$','$\hat{W}_{a,8}(t)$', '$\hat{W}_{a,9}(t)$',...
     '$\hat{W}_{a,10}(t)$','Interpreter','latex');
grid on

figure
p4 = plot(t, y(:,P.n+P.l+1:P.n+P.l+P.l));
title('Value Function Weights', 'FontSize',20,'Interpreter','latex')
xlabel('$Time (s)$','FontSize',20,'Interpreter','latex')
ylabel('$\hat{W}_{c}(t)$','FontSize',20,'Interpreter','latex')
legend('$\hat{W}_{,1}(t)$','$\hat{W}_{c,2}(t)$', '$\hat{W}_{c,3}(t)$',...
    '$\hat{W}_{c,4}(t)$','$\hat{W}_{c,5}(t)$', '$\hat{W}_{c,6}(t)$',...
    '$\hat{W}_{c,7}(t)$','$\hat{W}_{c,8}(t)$', '$\hat{W}_{c,9}(t)$',...
     '$\hat{W}_{c,10}(t)$','Interpreter','latex');;
grid on


function y_dot = closeLoopDynamics(t, y, P)

   s = y(1:P.n,1);
   x_0 = [-5;-5;5;5];
   s_0 =  [b(x_0(1),-7,5);b(x_0(2),-7,5);b(x_0(3),-5,7);b(x_0(4),-5,7)];
   theta_hat = y(P.n+P.l+P.l+P.l^2+1:P.n+P.l+P.l+P.l^2+P.q, 1);

   Y_s = reshape(y(P.n+P.l+P.l+P.l^2+P.q+1:P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q),:),P.n,P.q);
   Yf= reshape(y(P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+1:P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q),:),P.q,P.q);
   Gs = y(P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+1:P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+P.n,:);
   Xf =  y(P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+P.n+1:P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+P.n+P.q, :);
   x = y(P.n+P.l+P.l+P.l^2+P.q+(P.n*P.q)+(P.q*P.q)+P.n+P.q+1:end,:);
  

    
    bar_Y_f = 1000;
    update = 1;
    
    if norm(Yf) > bar_Y_f
        update = 0;
    end
    
     [u, Wa_hat_dot, Wc_hat_dot, Gamma_dot] = updateLawForActorCriticWeights(y, update, P);
     s_dot =  transFormedDynamics(t, s, u, P);
     x_dot = dynamics(t,x, u, P);
      Y_dot = Y(s,P)*update;  
      Yf_dot = (Y_s.'*Y_s)*update;
      Gs_dot = F1(s,P)+G(s,P)*u*update;
      Xf_dot = Y_s.'*(s-s_0-Gs)*update;
     
      gamma = 1000;
      theta_hat_dot = gamma*(Xf-Yf*theta_hat);

  
    y_dot = [s_dot;Wa_hat_dot;Wc_hat_dot;Gamma_dot(:);theta_hat_dot;...
        Y_dot(:);Yf_dot(:);Gs_dot;Xf_dot(:);x_dot];
    
end


function s_dot = transFormedDynamics(t, s, u, P)
   s_dot = F1(s,P)+F2(s,P)+G(s, P)*u;
end

function x_dot = dynamics(t, x, u, P)
     P.s2 = sin(x(2));
     P.c2 = cos(x(2));
     M = [P.p1+2*P.p3*P.c2 P.p2+P.p3*P.c2;P.p2+P.p3*P.c2 P.p2];
     V_m = [-P.p3*P.s2*x(4) -P.p3*P.s2*(x(3)+x(4));P.p3*P.s2*x(3) 0];
     D = diag([x(3);x(4);tanh(x(3));tanh(x(4))]);
     f_1x = [x(3);x(4);-M^-1*V_m*[x(3);x(4)]];
     f_2x = [0 0 0 0;0 0 0 0;-[M^-1 M^-1]*D];
     theta = [5.3;1.1;8.45;2.35];
     x_dot = f_1x+f_2x*theta+g(x,P)*u;
end 

function g_x = g(x,P)
    
     P.c2 = cos(x(2));
     M = [P.p1+2*P.p3*P.c2 P.p2+P.p3*P.c2;P.p2+P.p3*P.c2 P.p2];
     g_x = [0 0; 0 0; (M^-1).'];
end

function G_s = G(s,P)
    x = [binv(s(1),-7,5);binv(s(2),-7,5);binv(s(3),-5,7);binv(s(4),-5,7)];
    P.c2 = cos(x(2));
    M = [P.p1+2*P.p3*P.c2 P.p2+P.p3*P.c2;P.p2+P.p3*P.c2 P.p2];
    G_s = [0 0; 0 0; (M^-1).'];
end

function F2_s = F2(s,P)
    theta = [5.3;1.1;8.45;2.35];
    F2_s = Y(s, P)*theta;
end

function F1_s = F1(s,P)
     x = [binv(s(1),-7,5);binv(s(2),-7,5);binv(s(3),-5,7);binv(s(4),-5,7)];
     P.s2 = sin(x(2));
     P.c2 = cos(x(2));
     M = [P.p1+2*P.p3*P.c2 P.p2+P.p3*P.c2;P.p2+P.p3*P.c2 P.p2];
     V_m = [-P.p3*P.s2*x(4) -P.p3*P.s2*(x(3)+x(4));P.p3*P.s2*x(3) 0];
     F1_1 = B(s(1),-7,5)*x(3);
     F1_2 = B(s(2),-7,5)*x(4);
     F1_34 = (-M^-1*V_m)*[x(3);x(4)];
     F1_3 = B(s(3),-5,7)*F1_34(1);
     F1_4 = B(s(4),-5,7)*F1_34(2);
    F1_s = [F1_1;F1_2;F1_3;F1_4];
end

function y_s =  Y(s, P)

    x = [binv(s(1),-7,5);binv(s(2),-7,5);binv(s(3),-5,7);binv(s(4),-5,7)];
    P.c2 = cos(x(2));
    F_3_1 = B(s(3),-5,7)*(P.p2*x(3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_3_2 = B(s(3),-5,7)*-(x(4)*(P.p2+P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_3_3 = B(s(3),-5,7)*(P.p2*tanh(x(3)))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_3_4 = B(s(3),-5,7)*-(tanh(x(4))*(P.p2+P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_4_1 = B(s(4),-5,7)*-(x(3)*(P.p2 + P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_4_2 = B(s(4),-5,7)*(x(4)*(P.p1 + 2*P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_4_3 = B(s(4),-5,7)*-(tanh(x(3))*(P.p2+P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);
    F_4_4 = B(s(4),-5,7)*(tanh(x(4))*(P.p1+2*P.c2*P.p3))/(P.c2^2*P.p3^2+P.p2^2-P.p1*P.p2);

    y_s = [0 0 0 0; 0 0 0 0; F_3_1 F_3_2 F_3_3 F_3_4;F_4_1 F_4_2 F_4_3 F_4_4];
end

function [sigma_s, grad_sigma_s] = basis(s)
    s1 = s(1); s2=s(2); s3=s(3); s4=s(4);
    sigma_s = [s1*s3;s2*s4;s3*s2;s4*s1;s1*s2;s4*s3;s1^2;s2^2;s3^2;s4^2];
    grad_sigma_s = [[  s3,    0,   s1,    0]
                    [   0,   s4,    0,   s2]
                    [   0,   s3,   s2,    0]
                    [  s4,    0,    0,   s1]
                    [  s2,   s1,    0,    0]
                    [   0,    0,   s4,   s3]
                    [2*s1,    0,    0,    0]
                    [   0, 2*s2,    0,    0]
                    [   0,    0, 2*s3,    0]
                    [   0,    0,    0, 2*s4]];
end

function B_i = B(y, ai, Ai)
    B_i =(((ai^2*exp(y))-(2*ai*Ai)+(Ai^2*exp(-y)))/(Ai*ai^2-ai*Ai^2));
end

function b_y = b(y, ai, Ai)
     b_y  =  (log ((Ai/ai)*((ai-y)/(Ai-y))));  
end

function binv = binv(y, ai, Ai)
    binv =  ai*Ai*((exp(y))-1)/((ai*exp(y))-Ai);
end

function [u,Wa_hat_dot,Wc_hat_dot,Gamma_dot] = updateLawForActorCriticWeights(y, update, P)

   s = y(1:P.n, :); 
   Wa_hat = y(P.n+1:P.n+P.l, :); 
   Wc_hat = y(P.n+P.l+1:P.n+P.l+P.l, :);
   Gamma = reshape(y(P.n+P.l+P.l+1:P.n+P.l+P.l+P.l^2), P.l, P.l);
   theta_hat = y(P.n+P.l+P.l+P.l^2+1:P.n+P.l+P.l+P.l^2+P.q, 1); 
   
   summ_wc_r = zeros(size(Wc_hat));
   summ_gamma_r  = zeros(size(Wa_hat,1),size(Wa_hat,1));
   summ_wa_r  = zeros(size(Wa_hat));
    
   
    x1 = meshgrid(linspace(-4,3,10), linspace(-4,3,10));
    x2 = x1.'; x3 = x1; x4 = x1.';
    x_k = [x1(:) x2(:) x3(:) x4(:)];
    for k=1:P.N
        x_k_i = x_k(k, :);
        s_k = [b(x_k_i(1),-7,5) b(x_k_i(2),-7,5) b(x_k_i(3),-5,7)  b(x_k_i(4),-5,7)].';
        y_k = Y(s_k,P);
        G_k = G(s_k,P);
        F1_k = F1(s_k,P);
        [~, grad_sigma_sk] = basis(s_k);
        u_sk = uHat(s_k,Wa_hat,P);
        omega_k = grad_sigma_sk*(y_k*theta_hat+F1_k+G_k*u_sk);
        G_sigma_k = grad_sigma_sk*G_k*P.R^-1*G_k.'*grad_sigma_sk.';
        grad_sk_v_hat = Wc_hat.'*grad_sigma_sk;
        rho_k = 1+P.nu*(omega_k.'*omega_k);
        delta_hat_k = grad_sk_v_hat*(y_k*theta_hat+F1_k+G_k*u_sk)+(u_sk.'*P.R*u_sk)+(s_k.'*P.Q*s_k);
        summ_wc_r = summ_wc_r + (omega_k/rho_k)*delta_hat_k;
        summ_gamma_r = summ_gamma_r + (omega_k*omega_k.')/(rho_k.^2)*Gamma;
        summ_wa_r = summ_wa_r + P.kc2*G_sigma_k*Wa_hat*omega_k.'/(4*P.N*rho_k)*Wc_hat;
    end

   [~, grad_sigma_s] = basis(s);
   u =  uHat(s,Wa_hat,P);
   G_s = G(s,P);
   F1_s = F1(s,P);
   y_s = Y(s,P);
   omega = grad_sigma_s*(y_s*theta_hat+F1_s+G_s*u);
   G_sigma = grad_sigma_s*G_s*P.R^-1*G_s.'*grad_sigma_s.';
   grad_s_v_hat = Wc_hat.'*grad_sigma_s;
   rho = 1+P.nu*(omega.'*omega);
   delta_hat = grad_s_v_hat*(y_s*theta_hat+F1_s+G_s*u)+(u.'*P.R*u)+(s.'*P.Q*s);
   Wc_hat_dot = (-P.kc1*Gamma*(omega/rho)*delta_hat-P.kc2/P.N*Gamma*summ_wc_r)*update;
   
   Gamma_dot = (P.Beta*Gamma-(P.kc1*Gamma*((omega*omega.')/rho.^2)*Gamma)-...
       (P.kc2/P.N*Gamma*summ_gamma_r))*update;
   
   Wa_hat_dot = (-P.ka1*(Wa_hat-Wc_hat)-(P.ka2*Wa_hat)+((P.kc1*G_sigma.'...
       *Wa_hat*omega.')/(4*P.N*rho))*Wc_hat+summ_wa_r)*update;   
end


function u_hat = uHat(s, Wa_hat, P)
    [~, grad_sigma_s] = basis(s);
    u_hat = -1/2*P.R^-1*G(s,P).'*grad_sigma_s.'*Wa_hat;
end



