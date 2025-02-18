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
figure; %Fig 2(a)
plot(mu_range,Z1,LineWidth=2); hold on
plot(mu_range(idx), max_Z, 'r*', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(mu_range(idx)+0.02,max_Z,['($\mu_1^*,\ \phi_1(\mu_1^*)$)=(', num2str(mu_range(idx)),',',num2str(max_Z),')'],'Interpreter', 'latex')
xlabel('$\mu_1$')
ylabel('Fitness:$\phi_1(\mu_1)$')
fname = 'mu1phi1';
hfig = gcf;
%%
clear all; clc;close all
figure;%Fig 2(a)
a1=1.5;
a2=1;
b1=10;
b2=10;
c1=0.2;
c2=2;
fplot(@(x) a1.*sech(b1*x-c1).^2,'--',[0,1],LineWidth=2.5);hold on
fplot(@(x) a2.*sech(b2*x-c2).^2,':',[0,1],LineWidth=2.5);
xlabel('Death rates')
ylabel('Transmission function')
legend ('$\alpha$','$\beta$')
fname = 'alphaBeta';
hfig = gcf;
%%
%close all
syms mu mua
ep = 1/2.381; p = 0.956; omg = 1/3.119;
d = 0.000034;nu=0.031;%c=0.101;
% a1=1.5;
% a2=1;
% b1=25;
% b2=25;
% c1=2;
% c2=5;
% mua=0.5;
 alp=a1*((sech(b1*mua-c1)).^2);
% mua=d;
% alp=0.249;
A1=alp/(omg+mua);
A2=omg*(1-p)/(omg+mua);
beta=a2*((sech(b2*mu-c2)).^2);
R0=A1+A2*(beta)/(nu+mu);
% Convert the symbolic function R0 to a numeric function
f_numeric = matlabFunction(R0, 'Vars', [mu, mua]);


%%
close all
% Define the range for mu and fixed values for mua
mu_range = linspace(0, 1, 500);    % Range for mu
mua_values = [0.00001, 0.004008, 0.02]; % Different fixed values of mua

% Initialize the figure
figure;%Fig 3(c)
hold on;

% Define line styles for each iteration
line_styles = {'--', '-', ':', '-.'}; % Solid, dashed, dotted, dash-dot

% Loop over different fixed values of mua
for i = 1:length(mua_values)
    mua_fixed = mua_values(i); % Get the current fixed value of mua
    Z = f_numeric(mu_range, mua_fixed); % Compute Z for each fixed mua
    
    % Use the line style from the list based on the current iteration
    plot(mu_range, Z, 'DisplayName', ['$\mu_A$ = ', num2str(mua_fixed)], ...
        'LineWidth', 2, 'LineStyle', line_styles{i}); % Plot with line style
end
% Add labels, title, and legend
xlabel('$\mu_S$');
ylabel('Fitness: $\phi_2(\mu_A,\mu_S)$');
%title('Level Sets: \mu_a vs Z for Different \mu Values');
legend show;
legend(Location='best')
grid on;
hold off;
fname = 'muSVsPhi2';
hfig = gcf;
%%
close all
% Define the range for mua and fixed values for mu
mua_range = linspace(0, 1, 500);  % Range for mua
mu_values = [0.001, 0.17435, 0.5]; % Different fixed values of mu

% Initialize the figure
figure;%fig 3(b)
hold on;
% Define line styles for each iteration
line_styles = {'--', '-', ':', '-.'}; % Solid, dashed, dotted, dash-dot
% Loop over different fixed values of mu
for i = 1:length(mu_values)
    mu_fixed = mu_values(i); % Get the current fixed value of mu
    Z = f_numeric(mu_fixed, mua_range); % Compute Z for each fixed mu
    
    % Use the line style from the list based on the current iteration
    plot(mua_range, Z, 'DisplayName', ['$\mu_S$ = ', num2str(mu_fixed)], ...
        'LineWidth', 2, 'LineStyle', line_styles{i}); % Plot with line style
end

% Add labels, title, and legend
xlabel('$\mu_A$');
ylabel('Fitness: $\phi_2(\mu_A,\mu_S)$');
%title('Level Sets: \mu_a vs Z for Different \mu Values');
legend show;
grid on;
hold off;
fname = 'muAVsPhi2';
hfig = gcf;
%%
close all;
% Define ranges for mu and mua
mu_range = linspace(0, 1, 500);    % Range for mu
mua_range = linspace(0,1, 500);  % Range for mua

% Create meshgrid
[mu_mesh, mua_mesh] = meshgrid(mu_range, mua_range);

% Evaluate f_numeric on the meshgrid
Z = f_numeric(mu_mesh, mua_mesh);

% Find the maximum value of Z and its location
[max_Z, idx] = max(Z(:));
[max_mu_idx, max_mua_idx] = ind2sub(size(Z), idx);
max_mu = mu_mesh(max_mu_idx, max_mua_idx);
max_mua = mua_mesh(max_mu_idx, max_mua_idx);

% Plot the 3D surface plot
figure;%Fig 3(a)
surf(mu_mesh, mua_mesh, Z);
shading interp
hold on;
% Add a contour plot at the base of the surface
contour3(mu_mesh, mua_mesh, Z, 20, 'k') % '20' specifies the number of contour levels, 'k' for black lines
% Highlight the maximum value (argmax) with a red dot
plot3(max_mu, max_mua, max_Z, 'r*', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Add labels, title, and grid
xlabel('$\mu_S$');
ylabel('$\mu_A$');
zlabel('$\phi_2(\mu_A,\mu_S)$');
%title('3D Surface Plot of Z with Maximum Highlighted');
grid on;
hold off;

% Display the maximum value and location
disp(['Maximum Z occurs at \mu = ', num2str(max_mu), ' and \mu_a = ', num2str(max_mua)]);
disp(['Maximum Z value is ', num2str(max_Z)]);
fname = 'muAmuSVsPhi2';
hfig = gcf;
%%