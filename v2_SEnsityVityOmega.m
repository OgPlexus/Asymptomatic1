clear all; close all; clc;
%syms x y
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034;nu=0.031;
mu=0.17435;
%c=0.101;
a=1;
b=10;
c=2;
f = a * (sech(b*mu - c)).^2;       % Define f(x)
g = f / (mu + nu);
pri=linspace(1, 10, 100);
omg_values = 1./pri;
g1_max=g*ones(size(omg_values));
ess_mu=mu*ones(size(omg_values));
plot(1./omg_values,ess_mu,'r:','LineWidth', 2.5); hold on
%%
%close all; 
%clear all;clc
syms mu
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034;nu=0.031;%c=0.101;
a=1;
b=10;
c=2;

beta1=a*((sech(b*mu-c)).^2);
R0=(beta1)/(nu+mu);
% Convert the symbolic function R0 to a numeric function
f_numeric1 = matlabFunction(R0, 'Vars', mu);
mu_range = linspace(0, 1, 500);
Z1 = f_numeric1(mu_range);
[max_Z, idx] = max(Z1(:));
%clear all; close all; clc;
% Define the function f(x) and g(x)
syms mua omg
ep = 1/2.381; p = 0.956; %omg = 1/3.119;
d = 0.000034; %c = 0.101;
%nu=0.031;
mu=mu_range(idx);
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
% Create a range of y values
pri=linspace(1, 10, 100);
omg_values = 1./pri;%linspace(0, 1, 100); % 100 points in [0, 1]
mua_max_values = zeros(size(omg_values)); % Preallocate for solutions
g_max_values1 = zeros(size(omg_values)); % Preallocate for max g(x)

% Find argmax x for each y
for i = 1:length(omg_values)
    % Substitute the current y value into g(x)
    g_y = subs(g, omg, omg_values(i));
    
    % Maximize g(x) over the interval [0, 1]
    [mua_max, g_max] = fminbnd(@(mua_eval) -double(subs(g_y, mua, mua_eval)), 0, 1);
    
    % Store the results
    mua_max_values(i) = mua_max;       % x that maximizes g(x)
    g_max_values1(i) = -g_max;      % Maximum value of g(x)
end

% Plot y vs. x_max (argmax x)
%close all
plot(1./omg_values, mua_max_values, 'b--', 'LineWidth', 1.5);
xlabel('$\omega^{-1}$:  expected time in the asymptomatic state');
 ylabel('ESS: $\mu_S^*,\ \mu_A^*$');
% title('\mu_a^* vs. \nu');
legend('$\mu_S^*$','$\mu_A^*$','Location','best')
fname = 'omgVsVirulance';
hfig = gcf;

%%
 close all;
 figure;
plot(1./omg_values, mua_max_values, 'b', 'LineWidth', 1.5);
xlabel('$\omega^{-1}$:  expected time in the asymptomatic state');
ylabel('ESS: $\mu_A^*$');
grid on;
fname = 'OmgSenModel2MuA';
hfig = gcf;
%%
 close all;
 figure;
plot(1./omg_values, g1_max, 'r:', 'LineWidth', 2.5);hold on;
plot(1./omg_values, g_max_values1, 'b--', 'LineWidth', 1.5);
xlabel('$\omega^{-1}$:  expected time in the asymptomatic state');
ylabel('Maximum fitness value');
legend('Model 1: $\phi_1(\mu_1^*)$','Model 2: $\phi_2(\mu_A^*,\mu_S^*)$','Location','best')
grid on;
fname = 'OmgSenModelFit';
hfig = gcf;
