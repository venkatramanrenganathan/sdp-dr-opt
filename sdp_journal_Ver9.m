% Semidefinite Program to find probability bounds for one RV lying to the
% right of an hyperpane and another to its left half.
% Takes in two mean and variance data.
% Output: Sum of Probabilities = (False Positive + False Negative)
%
% problem data
clear all; close all; clc;
n = 2; % dimensions

x_position = 1:10;
optimum_values = zeros(length(x_position),1);
F_min = 2;
 
figure;

for i = 1 : length(x_position)

    % Legitimate Robot Position
    mu_leg = [i 5]'; % Distinct Neighbor

    % Malicious Robot Position
    mu_mal = [-i 5]';
    
    % Spoofed Robot Position
    mu_spoof = [-i 5]';

    % Covariance in X - corresponds to error in distance measurement
    % Covariance in Y - corresponds to error in angle measurement
    % Sigma = cov(randn(100,n));
    % Sigma = Sigma'*Sigma;

%     Sigma_leg = [1.5 0.5; 0.5 1];
%     Sigma_mal = [1.5 -0.5; -0.5 1];
    
    theta1 = atan2(mu_leg(2), mu_leg(1));
    Sigma_leg = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)]*[2 0; 0 1]*[cos(theta1) sin(theta1); -sin(theta1) cos(theta1)];
    theta2 = atan2(mu_mal(2), mu_mal(1));
    Sigma_mal = [cos(theta2) -sin(theta2); sin(theta2) cos(theta2)]*[2 0; 0 1]*[cos(theta2) sin(theta2); -sin(theta2) cos(theta2)];
    theta3 = theta2;
    Sigma_spoof = Sigma_mal;
    
    Sigma = blkdiag(Sigma_leg, Sigma_mal, Sigma_spoof);
    mu = [mu_leg;mu_mal;mu_spoof];

   
    A_1 = [-eye(n)  eye(n)   zeros(n)
           eye(n)   zeros(n) -eye(n)
           zeros(n) -eye(n)  eye(n)];
      
    A_2 = [-eye(n)  zeros(n) eye(n)
           zeros(n) eye(n)   -eye(n)
           eye(n)   -eye(n)  zeros(n)];
       
    A_3 = [eye(n)   -eye(n)  zeros(n)
           -eye(n)  eye(n)   zeros(n)
           zeros(n) zeros(n) zeros(n)]; 
       
    A_4 = [eye(n)   zeros(n) -eye(n)
           zeros(n) zeros(n) zeros(n)
           -eye(n)  zeros(n) eye(n)];
       
%% solve sdp
cvx_begin sdp
    variable Z(3*n,3*n) symmetric   % Legitimate  
    variable z(3*n,1)
    variable lambda
    maximize lambda
    subject to       
        trace(A_1*Z) >= 0;
        trace(A_2*Z) >= 0;
        trace(A_3*Z) <= F_min;
        trace(A_4*Z) <= F_min;
        [Z z; z' lambda] <= [Sigma + mu*mu' mu; mu' 1];
        [Z z; z' lambda] >= 0;      
cvx_end 

    optimum_values(i) = cvx_optval;

    h1 = plot_gaussian_ellipsoid_Ver1(1, i, mu_leg', Sigma_leg); % Format is (mu, sigma, sd) where sd = std deviation
    h2 = plot_gaussian_ellipsoid_Ver1(2, i, mu_mal', Sigma_mal); % sd = 1 by default.
    hold on;
    if (i == 1)
        plot(0,0,'x','color', 'k','Markersize', 10);
    end

end
 
%%

xlabel('x-position(cm)', 'interpreter', 'Latex');
ylabel('y-position(cm)', 'interpreter', 'Latex');
legend('Legitimate Neighbor', 'Malicious Neighbor','Receiving Robot');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
axis equal

figure;
distance = 2*x_position;
plot(distance, optimum_values, 'b');
xlabel('Distance between neighbors(cm)', 'interpreter', 'Latex');
ylabel('$$\lambda$$', 'interpreter', 'Latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);



