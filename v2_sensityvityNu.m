clear all; close all; clc;
% Define the function f(x) and g(x)
syms x y
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034;nu=0.031;%c=0.101;
a=1;
b=10;
c=2;
f = a * (sech(b*x - c)).^2;       % Define f(x)
g = f / (y + x);                % Define g(x)

% Create a range of y values
pri=linspace(1, 50, 100);
y_values = 1./pri;%linspace(0, 1, 100); % 100 points in [0, 1]
x_max_values = zeros(size(y_values)); % Preallocate for solutions
g_max_values1 = zeros(size(y_values)); % Preallocate for max g(x)

% Find argmax x for each y
for i = 1:length(y_values)
    % Substitute the current y value into g(x)
    g_y = subs(g, y, y_values(i));
    
    % Maximize g(x) over the interval [0, 1]
    [x_max, g_max] = fminbnd(@(x_eval) -double(subs(g_y, x, x_eval)), 0, 1);
    
    % Store the results
    x_max_values(i) = x_max;       % x that maximizes g(x)
    g_max_values1(i) = -g_max;      % Maximum value of g(x)
end

% Plot y vs. x_max (argmax x)
%close all
figure;
plot(1./y_values, x_max_values, 'k', 'LineWidth', 1.5);
xlabel('$\nu^{-1}$: mean period in symptomatic stage');
ylabel('ESS: $\mu_1^*$');
%title('x = argmax g(x) vs. y');
grid on;
fname = 'senNuModel1';
hfig = gcf;
hold on
%%
syms mu mua nu

% Define constants
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034; %c = 0.101;
% a1 = 2; a2 = 2;
% b1 = 5; b2 = 5;
a1=1.5;
a2=1;
b1=10;
b2=10;
c1=0.2;
c2=2;

% Define functions
alp = a1 * (sech(b1 * mua - c1))^2;
A1 =  alp /  (omg + mua);
A2 = omg * (1 - p) / (omg + mua);
beta = a2 * (sech(b2 * mu - c2))^2;
g = A1 + A2 * (beta) / (nu + mu);

% Create a range of nu values
pri=linspace(1, 50, 100);
nu_values = 1./pri;%linspace(0, 1, 100); % 100 points in [0, 1]
%nu_values = linspace(0, 1, 100);%linspace(0, 1, 100)
mu_max_values = zeros(size(nu_values));  % Preallocate for mu*
mua_max_values = zeros(size(nu_values)); % Preallocate for mua*
g_max_values = zeros(size(nu_values));   % Preallocate for max g

% Find argmax (mu*, mua*) for each nu
for i = 1:length(nu_values)
    % Substitute the current nu value into g(mu, mua)
    g_nu = subs(g, nu, nu_values(i));

    % Maximize g(mu, mua) over mu and mua
    g_handle = matlabFunction(g_nu, 'Vars', [mu, mua]); % Convert to function handle
    [x_opt, g_max] = fmincon(@(x) -g_handle(x(1), x(2)), ... % Negative for maximization
                              [0.5, 0.5], ... % Initial guess
                              [], [], [], [], ... % No linear constraints
                              [0, 0], [1, 1]); % Bounds for mu and mua

    % Store the results
    mu_max_values(i) = x_opt(1);       % mu*
    mua_max_values(i) = x_opt(2);      % mua*
    g_max_values(i) = -g_max;          % Maximum value of g
end

% Plot nu vs. mu*
% figure;
 plot(1./nu_values, mu_max_values, 'r.', 'LineWidth', 1.5);
% xlabel('\nu');
% ylabel('\mu^*');
% %title('\mu^* vs. \nu');
% grid on;

% Plot nu vs. mua*
% figure;
plot(1./nu_values, mua_max_values, 'b--', 'LineWidth', 1.5);
% xlabel('\nu');
% ylabel('\mu_a^*');
% title('\mu_a^* vs. \nu');
legend('$\mu_1^*$','$\mu_S^*$','$\mu_A^*$','Location','best')
fname = 'nuVsVirulance';
hfig = gcf;
%%max value
 % plot(nu_values, g_max_values, 'r-', 'LineWidth', 1.5);

 %%
 close all;
 figure;
plot(1./nu_values, mua_max_values, 'b', 'LineWidth', 1.5);
xlabel('$\nu^{-1}$: mean period in symptomatic stage');
ylabel('ESS: $\mu_A^*$');
grid on;
fname = 'senNuModel2MuA';
hfig = gcf;

 %%
 close all;
 figure;
plot(1./nu_values, mu_max_values, 'r', 'LineWidth', 1.5);
xlabel('$\nu^{-1}$: mean period in symptomatic stage');
ylabel('ESS: $\mu_S^*$');
grid on;
fname = 'senNuModel2MuS';
hfig = gcf;
%%
 close all;
 figure;
plot(1./nu_values, g_max_values1, 'r:', 'LineWidth', 2);hold on;
plot(1./nu_values, g_max_values, 'b--', 'LineWidth', 1.5);
xlabel('$\nu^{-1}$: mean period in symptomatic stage');
ylabel('Maximum fitness value');
legend('Model 1: $\phi_1(\mu_1^*)$','Model 2: $\phi_2(\mu_A^*,\mu_S^*)$','Location','best')
grid on;
fname = 'senNuModelFit';
hfig = gcf;