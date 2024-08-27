% 
% This script is designed to tune the following six control parameters.
%   k1, k2, kappa_bar, k3s, k4, D
% 
% email: wusixstx@163.com
% date:  2024.08.18


clear all; clc; close all;


% Given a specific parameter k4 
k4         = 30;
kappa_bar  = tand(17);

% Calculate temporary parameters, which denpend on the WIP model.
[G_ubar, G_bar] = get_G(k4, kappa_bar);
k_alpha = get_kalpha();

% Condition required by the tracking controller (see Proposition 2)
k3 = k4^2 / 4;
k2 = k4 / 25 / G_bar;
k_ztilde2 = 3/(k2*G_ubar);

% Condition required by the Theorem (see the small gain condition in its proof)
c = 1.01;
con(1) = k_ztilde2*(3*(k_alpha+1)*c);
con(2) = k_ztilde2*(3*(k_alpha+1)*c)*k_alpha*k4*k2/2;
con(3) = k_ztilde2*(3*(k_alpha+1)^2*c^2);
con(4) = k_ztilde2*(3*(k_alpha+1)^2*c^2)*k2*k4/2;

k      = max(con)*c; % k>max(con);
k1     = (k_alpha + 1) * c / k;

% Print result
fprintf('k1: %f, (alpha(s)=k1*s)\n', k1);
fprintf('k2: %f\n', k2);
fprintf('k3: %f\n', k3);
fprintf('k4: %f\n', k4);
fprintf('kappa_bar: %f\n', kappa_bar);




function [G_ubar, G_bar, z3z4_ubar, z3z4_bar] = get_G(k4, kappa_bar)

    % Define the objective functions with penalty terms
    G_new_ubar = @(x) G(x(1), x(2)) + (abs(x(1)) > 2 * kappa_bar) * 1e10 + (abs(x(2)) > k4 * kappa_bar / 2) * 1e10;
    G_new_bar  = @(x) -G(x(1), x(2)) + (abs(x(1)) > 2 * kappa_bar) * 1e10 + (abs(x(2)) > k4 * kappa_bar / 2) * 1e10;

    % Find the minimum values
    [z3z4_ubar, G_ubar] = fminsearch(G_new_ubar, [0; 0]);
    [z3z4_bar, G_bar]   = fminsearch(G_new_bar, [0; 0]);
    G_bar = -G_bar;  % Correct the sign of G_bar

end



function k_alpha = get_kalpha()
phi_zw_minus = @(x) (-abs(phizw(x(1), x(2))));
[z1z3, k_alpha_minus] = fminsearch(phi_zw_minus, [0;0]);
k_alpha = -k_alpha_minus;
end

function dphizw_dz3 = phizw(z1, z3)
m_w = 5.00;                     % Mass of the wheel (kg)
r_w = 0.08;                     % Radius of the wheel (m)
m_b = 60.0;                     % Mass of the body (kg)
l_b = 0.40;                     % Half-length of the rod (m)

g   = 9.8;                      % Acceleration due to gravity (m/s^2)
L_b = l_b * 2.0;                % Total length of the rod (m)
J_w = 0.5 * m_w * r_w^2;        % Moment of inertia of the wheel (cylindrical), 1/2*m*r^2
J_b = 1/3 * m_b * l_b^2;        % Moment of inertia of the rod (thin rod), 1/3*m*l^2

C1 = J_w + (m_b + m_w) * r_w^2;
C2 = m_b *  r_w * l_b;
C3 = J_b +  m_b * l_b^2;

y      = z1 - psi(atan(z3)) + atan(z3);
% \partial phizw/ \partial z3 = (-\partial psi / \partial atan(z3) + 1) \partial atan(z3)/ \partial z3;
dphizw_dz3 = (- psi_d(atan(z3))+1 ) * 1/(1+z3^2);

    function res = psi(qs)
        res = 2*qs - log( (C2 + C1*cos(qs) + sin(qs)*sqrt(C2^2 - C1^2))...
                            /(C1 + C2*cos(qs)) ) ...
                            * (C1 - C3) / sqrt(C2^2 - C1^2);
    end 
    
    function res = psi_d(qs)
        res = 2 - (C1 - C3)/(C1 + C2*cos(qs)); 
    end

end

function y = G(z3, z4)
m_w = 5.00;                     % Mass of the wheel (kg)
r_w = 0.08;                     % Radius of the wheel (m)
m_b = 60.0;                     % Mass of the body (kg)
l_b = 0.40;                     % Half-length of the rod (m)

g   = 9.8;                      % Acceleration due to gravity (m/s^2)
L_b = l_b * 2.0;                % Total length of the rod (m)
J_w = 0.5 * m_w * r_w^2;        % Moment of inertia of the wheel (cylindrical), 1/2*m*r^2
J_b = 1/3 * m_b * l_b^2;        % Moment of inertia of the rod (thin rod), 1/3*m*l^2

C1 = J_w + (m_b + m_w) * r_w^2;
C2 = m_b *  r_w * l_b;
C3 = J_b +  m_b * l_b^2;

y =  C2 * g / r_w / (C1*sqrt(1+z3^2)+C2) ...
        + (C2*(C2+C3*sqrt(z3^2+1))*z4^2) ...
        / (1+z3^2)^2 / (C2+C1*sqrt(1+z3^2))^2;
end