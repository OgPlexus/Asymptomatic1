close all; clear all;clc
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
mu_range(idx)
plot(mu_range,Z1,LineWidth=2); hold on
plot(mu_range(idx), max_Z, 'r*', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(mu_range(idx)+0.02,max_Z,['($\mu_1^*,\ \phi_1(\mu_1^*)$)=(', num2str(mu_range(idx)),',',num2str(max_Z),')'],'Interpreter', 'latex')
xlabel('$\mu_1$')
ylabel('Fitness:$\phi_1(\mu_1)$')
fname = 'mu1phi1';
hfig = gcf;
%%
syms mu mua
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034;nu=0.031;%c=0.101;
a1=1.5;
a2=1;
b1=10;
b2=10;
c1=0.2;
c2=2;
 alp=a1*((sech(b1*mua-c1)).^2);
% mua=d;
% alp=0.249;
A1=ep*alp/((ep+d)*(omg+mua));
A2=ep*omg*(1-p)/((ep+d)*(omg+mua));
beta=a2*((sech(b2*mu-c2)).^2);
R0=A1+A2*(beta)/(nu+mu);
mu_range = linspace(0, 1, 500);
% Convert the symbolic function R0 to a numeric function
f_numeric = matlabFunction(R0, 'Vars', [mu, mua]);
close all;
figure;
plot(mu_range,Z1,LineWidth=2); hold on
plot(mu_range, f_numeric (mu_range,0))
xlabel('$\mu_1$ or $\mu_S$')
ylabel('Fitness:$\phi_1(\mu_1)$ or $\phi_2(0, \mu_s)$')

%%
fplot(R0,[0,1],LineWidth=2)
%hold on
%%
legend('$\mu=0$','$\mu=0.25$','$\mu=0.5$','$\mu=1$',Location='best')

%%
xlabel('Regression Coefficient: $a$')
ylabel('Basic reproduction number $\mathcal{R}_0$')
fname = 'rmuVsR0';
hfig = gcf;
%%
clear all; close all;
syms a mu
b=0.001;d=0.001;
c=0.01;
R0=b*(a*mu+c)/(d*(d+mu));
loweBound=-c
%
f_numeric = matlabFunction(R0, 'Vars', [mu, a]);
[X, Y] = meshgrid(linspace(0.1, 1, 50), linspace(-1, 1, 50));
Z = f_numeric(X, Y);
% Plot the surface
surf(X, Y, Z)

% Label the axes
xlabel('mu')
ylabel('r')
zlabel('R0(\mu, r)')

% Add shading for better visualization
shading interp

%%
clear all; close all
% Define symbolic variables
syms a mu

% Define constants
nu=0.1;
d=0.001;
c=0.101;

% Define the symbolic function for R0

R0=(a*mu+c)/(d+nu+mu);

% Convert the symbolic function R0 to a numeric function
f_numeric = matlabFunction(R0, 'Vars', [mu, a]);

% Define mu values ranging from 0.0001 to 1
mu_vals = linspace(0, 1, 100);

% Initialize the grid for 'a' and 'Z' (R0 values)
a_grid = [];
Z = [];

% Loop over each mu value
for i = 1:length(mu_vals)
    mu_val = mu_vals(i);
    
    % Compute the lower bound for 'a' for the given 'mu'
    lowerBound = max(-c / mu_val,-2);
    
    % Define 'a' values for the current 'mu' value (from lowerBound to 1)
    a_vals = linspace(lowerBound, 2, 100);
    
    % Store the 'a' values (different for each 'mu')
    a_grid = [a_grid; a_vals];
    
    % Evaluate the function for the current 'mu' and corresponding 'a' values
    Z_row = f_numeric(mu_val * ones(size(a_vals)), a_vals);
    
    % Store the computed 'R0' values in the Z matrix
    Z = [Z; Z_row];
end

% Create a meshgrid for mu and a to use in pcolor
[mu_mesh, a_mesh] = meshgrid(mu_vals, linspace(min(a_grid(:)), 2, size(Z, 1)));

% Plot the pcolor graph
surf(mu_mesh, a_mesh, Z)
shading interp

% Label the axes
xlabel('$\mu$')
ylabel('$a$')
zlabel('$\mathcal{R}_0$')
%title('Pcolor plot of R_0 as a function of \mu and a')

% Add colorbar for reference
%colorbar
%%
clear all; clc
syms a 
alp=0.249;ep=0.42;p=0.956;omg=0.3206;
c=0.101;d=0.001;nu=0.1;
A1=ep*alp/((ep+d)*(omg+d));
A2=ep*omg*(1-p)/((ep+d)*(omg+d));
beta=a*mu+c;
mu=1
R0=A1+A2*(beta)/(d+nu+mu);
%%
fplot(R0,[max(-c/mu,-2),3],LineWidth=2)
hold on

%%
clear all; clc; close all
syms a mu
alp=0.249;ep=0.42;p=0.956;omg=0.3206;
c=0.101;d=0.001;nu=0.1;
A1=ep*alp/((ep+d)*(omg+d));
A2=ep*omg*(1-p)/((ep+d)*(omg+d));
R0=A1+A2*(a*mu+c)/(d+nu+mu);
% Convert the symbolic function R0 to a numeric function
f_numeric = matlabFunction(R0, 'Vars', [mu, a]);

% Define mu values ranging from 0.0001 to 1
mu_vals = linspace(0, 1, 100);

% Initialize the grid for 'a' and 'Z' (R0 values)
a_grid = [];
Z = [];

% Loop over each mu value
for i = 1:length(mu_vals)
    mu_val = mu_vals(i);
    
    % Compute the lower bound for 'a' for the given 'mu'
    lowerBound = max(-c / mu_val,-2);
    
    % Define 'a' values for the current 'mu' value (from lowerBound to 1)
    a_vals = linspace(lowerBound, 2, 100);
    
    % Store the 'a' values (different for each 'mu')
    a_grid = [a_grid; a_vals];
    
    % Evaluate the function for the current 'mu' and corresponding 'a' values
    Z_row = f_numeric(mu_val * ones(size(a_vals)), a_vals);
    
    % Store the computed 'R0' values in the Z matrix
    Z = [Z; Z_row];
end

% Create a meshgrid for mu and a to use in pcolor
[mu_mesh, a_mesh] = meshgrid(mu_vals, linspace(min(a_grid(:)), 2, size(Z, 1)));

% Plot the pcolor graph
surf(mu_mesh, a_mesh, Z)
shading interp

% Label the axes
xlabel('$\mu$')
ylabel('$a$')
zlabel('$\mathcal{R}_0$')
%title('Pcolor plot of R_0 as a function of \mu and a')

% Add colorbar for reference
%%
syms mu
nu=0.1;d=0.001;
c=0.101;a=2
R0=(a*mu+c)/(d+nu+mu);
%%
fplot(R0,[0,1],LineWidth=2)
hold on
%%
legend('$a=-0.2$','$a=0.2$','$a=2$',Location='best')
xlabel('Virulence: $\mu$')
ylabel('Basic reproduction number $\mathcal{R}_0$')
fname = 'muVsR0';
hfig = gcf;
%%
%%
clear all; clc
syms mu mua
alp=0.249;ep=0.42;p=0.956;omg=0.3206;
c=0.101;d=0.001;nu=0.1;
a=2
% mua=0.5;
 alp=5*sech(25*mua-2);
% mua=d;
% alp=0.249;
A1=ep*alp/((ep+d)*(omg+mua));
A2=ep*omg*(1-p)/((ep+d)*(omg+mua));
beta=a*sech(25*mu-2);
R0=A1+A2*(beta)/(d+nu+mu);
%%
fplot(R0,[0,1],LineWidth=2)
hold on
%%
legend('$a=-0.3$','$a=0.8$','$a=2$',Location='best')
xlabel('Virulence: $\mu$')
ylabel('Basic reproduction number $\mathcal{R}_0$')
fname = 'muVsR0Covid';
hfig = gcf;