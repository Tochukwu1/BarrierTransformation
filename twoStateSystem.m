close all; clear all; clc;

% Simulation Parameters
P.Q = 10*eye(2);
P.R = 0.1;
P.N = 100;
P.l = 3; 
P.n = 2;
P.kc = 5;
P.ka1 = 1500;
P.ka2 = 0.01;
P.Beta = 1;
P.nu = 1;
P.Beta_1 = 5;
P.k = 10;
P.alpha = 1;


% Initial Conditions
x_0 = [-6;6];
x_hat_0 = [-6;4];
s_0 = [b(-6,-7,5) b(6,-5,7)].';
s_hat_0 = [b(-6,-7,5) b(4,-5,7)].';
gamma_0 = eye(3);
wc_hat0 = [10;1/2;1/2];
wa_hat0 =  wc_hat0;

y0 = [s_0;s_hat_0;wa_hat0;wc_hat0;gamma_0(:);x_0;x_hat_0;0];
tspan = [0 8];

% Simulation

options = odeset('OutputFcn',@odeplot, 'OutputSel', 1:P.n+P.n);
[t,y] = ode45(@(t,y) closeLoopDynamics(t,y, P), tspan, y0, options);

% Plots
x = y(:, P.n+P.n+P.l+P.l+P.l^2+1:P.n+P.n+P.l+P.l+P.l^2+P.n);
x_hat = y(:, P.n+P.n+P.l+P.l+P.l^2+P.n+1:P.n+P.n+P.l+P.l+P.l^2+P.n+P.n);
x_tilde = x - x_hat;
plot(t, x_tilde);
legend('$\tilde{x_1}$', '$\tilde{x_2}$', 'Interpreter','Latex')
title('Graph of $\tilde{x}$', 'Interpreter','Latex')
grid on

figure
s = y(:, 1:P.n);
plot(s(:,1), s(:,2));
title('Graph of s', 'Interpreter','Latex')
grid on

figure
plot(t, y(:, P.n+P.n+1:P.n+P.n+P.l));
title('Policy Weights', 'FontSize',20,'Interpreter','latex')
xlabel('$Time (s)$','FontSize',20,'Interpreter','latex')
ylabel('$\hat{W}_{a}(t)$','FontSize',20,'Interpreter','latex')
legend('$\hat{W}_{a,1}(t)$','$\hat{W}_{a,2}(t)$', '$\hat{W}_{a,3}(t)$', 'Interpreter','latex');
grid on

figure
plot(t, y(:,P.n+P.n+P.l+1:P.n+P.n+P.l+P.l));
title('Value Function Weights', 'FontSize',20,'Interpreter','latex')
xlabel('$Time (s)$','FontSize',20,'Interpreter','latex')
ylabel('$\hat{W}_{c}(t)$','FontSize',20,'Interpreter','latex')
legend('$\hat{W}_{c,1}(t)$','$\hat{W}_{c,2}(t)$', '$\hat{W}_{c,3}(t)$', 'Interpreter','latex');
grid on


function y_dot = closeLoopDynamics(t, y, P)

    s = y(1:P.n, :);
    s_hat = y(P.n+1:P.n+P.n, :);
    x = y(P.n+P.n+P.l+P.l+P.l^2+1:P.n+P.n+P.l+P.l+P.l^2+P.n, :);
    x_hat = y(P.n+P.n+P.l+P.l+P.l^2+P.n+1:P.n+P.n+P.l+P.l+P.l^2+P.n+P.n, :);

    nu_1 = feedbackTerm(y, P);
    [u ,Wa_hat_dot,Wc_hat_dot,Gamma_dot]= updateLawForActorCriticWeights(y,P);

    eta_dot = dynamicFilterOutput(t,y, u, nu_1, P);
    x_dot = dynamics(t,x,u);
    x_hat_dot = dynamicsEstimate(t, x_hat, u, nu_1);
    s_hat_dot = transFormedDynamicsEstimate(t, s_hat, u, nu_1);
    s_dot = transFormedDynamics(t, s, u);

    y_dot = [s_dot;s_hat_dot;Wa_hat_dot;Wc_hat_dot;Gamma_dot(:);x_dot;x_hat_dot;eta_dot];
end

function s_dot = transFormedDynamics(t, s, u) 
   s1_dot = H(s);
   s2_dot = F(s)+G(s)*u;
   s_dot = [s1_dot;s2_dot];
end

function x_dot = dynamics(t, x, u)
   x1_dot = x(2);
   x2_dot = f(x)+g(x)*u;
   x_dot = [x1_dot;x2_dot];
end 

function x_hat_dot = dynamicsEstimate(t, x_hat, u, nu_1)
   x1_hat_dot = x_hat(2);
   x2_hat_dot = f(x_hat)+g(x_hat)*u + nu_1;
   x_hat_dot = [x1_hat_dot;x2_hat_dot];
end 

function f_x = f(x)
    f_x = -x(1)-1/2*x(2)*(1-(cos(2*x(1))+2)^2);
end 

function g_x = g(x)
     g_x = cos(2*x(1))+2;
end

function H_s = H(s)
    H_s = B(s(1),-7,5)*binv(s(2),-5,7);
end

function G_s = G(s)
    x = [binv(s(1),-7,5); binv(s(2),-5,7)];
    G_s = B(s(2),-5,7)*(cos(2*x(1))+2);
end

function F_s = F(s)
    x = [binv(s(1),-7,5); binv(s(2),-5,7)];
    F_s = B(s(2),-5,7)*(-x(1)-1/2*x(2)*(1-(cos(2*x(1))+2)^2));
end

function [sigma_s_hat, grad_sigma_s_hat] = basis(s_hat)

    sigma_s_hat = [s_hat(1)^2;s_hat(1)*s_hat(2);s_hat(2)^2];
    grad_sigma_s_hat = [2*s_hat(1) 0;s_hat(2) s_hat(1);0 2*s_hat(2)];
end

function B_i = B(y, ai, Ai)
    B_i =(((ai^2*exp(y))-(2*ai*Ai)+(Ai^2*exp(-y)))/(Ai*ai^2-ai*Ai^2));
end

function b_y = b(y, ai, Ai)
     b_y  =  log((Ai/ai)*((ai-y)/(Ai-y)));  
end

function binv = binv(y, ai, Ai)
    binv =  ai*Ai*((exp(y))-1)/((ai*exp(y))-Ai);
end

function nu_1 = feedbackTerm(y, P)
    x = y(P.n+P.n+P.l+P.l+P.l^2+1:P.n+P.n+P.l+P.l+P.l^2+P.n, :);
    x_hat = y(P.n+P.n+P.l+P.l+P.l^2+P.n+1:P.n+P.n+P.l+P.l+P.l^2+P.n+P.n, :);
    eta = y(P.n+P.n+P.l+P.l+P.l^2+P.n+P.n+1:end,:);
    nu_1 = (P.alpha^2*(b(x(1),-7,5)-b(x_hat(1),-7,5))-(P.k+P.alpha+P.Beta_1)*eta)/...
        B(b(x_hat(1),-7,5), -7,5);
end

function eta_dot = dynamicFilterOutput(t,y, u, nu_1, P)
    s = y(1:P.n, :);
    s_hat = y(P.n+1:P.n+P.n, :);
    eta = y(P.n+P.n+P.l+P.l+P.l^2+P.n+P.n+1:end,:);
    
    s_dot = transFormedDynamics(t, s, u);
    s_hat_dot = transFormedDynamicsEstimate(t, s_hat, u, nu_1);
    r_j = s_dot(1)-s_hat_dot(1) +P.alpha*(s(1)-s_hat(1))+eta;
    eta_dot = -P.Beta_1*eta-P.k*r_j-P.alpha*(s_dot(1)-s_hat_dot(1));
end

function s_hat_dot = transFormedDynamicsEstimate(t, s_hat, u, nu_1)
      
      nu_2 = B(s_hat(2),-5,7)*nu_1;
      s1_hat_dot = H(s_hat);
      s2_hat_dot = F(s_hat)+G(s_hat)*u + nu_2;
      s_hat_dot = [s1_hat_dot;s2_hat_dot];
end

function u_hat = uHat(s_hat, Wa_hat, P)
    [~, grad_sigma_s_hat] = basis(s_hat);

    u_hat = -1/2*P.R^-1*G(s_hat).'*grad_sigma_s_hat(:,2).'*Wa_hat;
end


function [u_hat,Wa_hat_dot,Wc_hat_dot,Gamma_dot] = updateLawForActorCriticWeights(y, P)

   s_hat = y(P.n+1:P.n+P.n, :);
   Wa_hat = y(P.n+P.n+1:P.n+P.n+P.l, :);
   Wc_hat = y(P.n+P.n+P.l+1:P.n+P.n+P.l+P.l, :);
   Gamma = reshape(y(P.n+P.n+P.l+P.l+1:P.n+P.n+P.l+P.l+P.l^2), P.l, P.l);
   
   summ_wc_r = zeros(size(Wc_hat));
   summ_gamma_r  = zeros(size(Wa_hat,1),size(Wa_hat,1));
   summ_wa_r  = zeros(size(Wa_hat));
    
   
    x1_hat = meshgrid(linspace(-2,2,10), linspace(-2,2,10));
    x2_hat = x1_hat.';
    x_hat_k = [x1_hat(:) x2_hat(:)];
    for k=1:P.N
        x_k_hat_i = x_hat_k(k, :);
        s_hat_k = [b(x_k_hat_i(1),-7,5) b(x_k_hat_i(2),-5,7)].';
        G_k = G(s_hat_k);
        [~, grad_sigma_hat_sk] = basis(s_hat_k);
        u_hat_sk = uHat(s_hat_k,Wa_hat,P);
        H_s_hat_k = H(s_hat_k);
        F_k = F(s_hat_k);
        Qsk_prime = Q(s_hat_k, P);
        omega_k=grad_sigma_hat_sk(:,1)*H_s_hat_k+grad_sigma_hat_sk(:,2)*(F_k+G_k*u_hat_sk);
        G_sigma_k = grad_sigma_hat_sk(:, 2)*G_k*P.R^-1*G_k.'*grad_sigma_hat_sk(:, 2).';
        grad_sk_v_hat = Wc_hat.'*grad_sigma_hat_sk;
        rho_k = 1+P.nu*(omega_k.'*omega_k);
        delta_hat_k = grad_sk_v_hat(:,1)*H_s_hat_k+Qsk_prime+grad_sk_v_hat(:,2)*...
            (F_k+G_k*u_hat_sk)+(u_hat_sk.'*P.R*u_hat_sk);
        summ_wc_r = summ_wc_r + (omega_k/rho_k)*delta_hat_k;
        summ_gamma_r = summ_gamma_r + (omega_k*omega_k.')/(rho_k^2)*Gamma;
        summ_wa_r = summ_wa_r + P.kc*G_sigma_k*Wa_hat*omega_k.'/(4*P.N*rho_k)*Wc_hat;
    end

   u_hat =  uHat(s_hat,Wa_hat,P);

   Wc_hat_dot = -P.kc/P.N*Gamma*summ_wc_r;
   
   Gamma_dot = P.Beta*Gamma-(P.kc/P.N*Gamma*summ_gamma_r);
   
   Wa_hat_dot = -P.ka1*(Wa_hat-Wc_hat)-(P.ka2*Wa_hat)+summ_wa_r;   
end

function Qs_prime = Q(s, P)
    Qs_prime = s.'*P.Q*s;
end