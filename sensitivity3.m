%clear all;
syms mu nu

% Define constants
ep = 0.42; p = 0.956; omg = 0.3206;
c = 0.101; d = 0.001;mua=0.001;
% a1 = 10; a2 = 2;
% b1 = 5; b2 = 5;

% Define functions
alp = a1 * (sech(b1 * mu - 2))^2;
A1 = ep * alp / ((ep + d) * (omg + mua));
A2 = ep * omg * (1 - p) / ((ep + d) * (omg + mua));
beta = a2 * (sech(b2 * mu - 2))^2;
g = A1 + A2 * (beta) / (nu + mu);

% Create a range of y values
y_values = linspace(0, 1, 100); % 100 points in [0, 1]
x_max_values = zeros(size(y_values)); % Preallocate for solutions
g_max_values = zeros(size(y_values)); % Preallocate for max g(x)

% Find argmax x for each y
for i = 1:length(y_values)
    % Substitute the current y value into g(x)
    g_y = subs(g, nu, y_values(i));
    
    % Maximize g(x) over the interval [0, 1]
    [x_max, g_max] = fminbnd(@(x_eval) -double(subs(g_y, mu, x_eval)), 0, 1);
    
    % Store the results
    x_max_values(i) = x_max;       % x that maximizes g(x)
    g_max_values(i) = -g_max;      % Maximum value of g(x)
end

% Plot nu vs. mu*
% figure;
 plot(y_values, x_max_values, 'o-', 'LineWidth', 1.5);
% xlabel('\nu');
% ylabel('\mu^*');
% %title('\mu^* vs. \nu');
% grid on;

% Plot nu vs. mua*
% figure;
%plot(nu_values, mua_max_values, 'b-', 'LineWidth', 1.5);
% xlabel('\nu');
% ylabel('\mu_a^*');
% title('\mu_a^* vs. \nu');
% grid on;

%%max value
 % plot(nu_values, g_max_values, 'r-', 'LineWidth', 1.5);