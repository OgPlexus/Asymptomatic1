clear all; clc;close all
figure;
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
% Define ranges for mu and mua
mu_range = linspace(0, 1, 500);    % Range for mu
mua_range = linspace(0,1, 500);  % Range for mua

% Create meshgrid
[mu_mesh, mua_mesh] = meshgrid(mu_range, mua_range);

% Evaluate f_numeric on the meshgrid
Z = f_numeric(mu_mesh, mua_mesh);
%%
% Plot the surface
surf(mu_mesh, mua_mesh, Z)
shading interp
hold on

% Add a contour plot at the base of the surface
contour3(mu_mesh, mua_mesh, Z, 20, 'k') % '20' specifies the number of contour levels, 'k' for black lines
% Label the axes
xlabel('$\mu_S$')
ylabel('$\mu_A$')
zlabel('Fitness: $\phi_2$')
%title('Pcolor plot of R_0 as a function of \mu and a')
%%
% Find the maximum value in Z and its index
[maxValue, linearIndex] = max(Z(:));

% Convert the linear index to subscripts (row, col)
[row, col] = ind2sub(size(Z), linearIndex);

% Get the corresponding mu and mua values
mu_max = mu_mesh(row, col);
mua_max = mua_mesh(row, col);

% Display the results
disp(['Maximum value of f_numeric: ', num2str(maxValue)])
disp(['Argmax coordinates: mu = ', num2str(mu_max), ', mua = ', num2str(mua_max)])
%%
contour(mu_mesh,mua_mesh,Z);
shading interp;
colorbar;
%clim([0,3])
%%
% Define ranges for mu and specific values for mua
mu_range = linspace(0, 1, 100);  % Range for mu
mua_values = [ 0.001, 0.06,0.08,0.1];     % Specific mua values

% Preallocate for results
f_numeric_values = zeros(length(mu_range), length(mua_values));

% Loop through mua values and compute f_numeric
for i = 1:length(mua_values)
    mua_fixed = mua_values(i);  % Current mua value
    f_numeric_values(:, i) = f_numeric(mu_range, mua_fixed); % Compute for fixed mua
end

% Plot the results
close all
figure;
hold on;
for i = 1:length(mua_values)
    plot(mu_range, f_numeric_values(:, i), 'DisplayName', ['$\mu_A$ = ', num2str(mua_values(i))],'LineWidth',2);
end
hold off;

% Add labels, legend, and title
xlabel('$\mu$')
ylabel('fitness:$\phi_2$')
%title('2D Plots of \mu vs f_{numeric} at Different \mu_a Values')
legend show
legend('Location','best')
grid on

%%
% Define range for mua and specific values for mu
mua_range = linspace(0.01, 1, 5000); % Range for mua
mu_values = [0.01, 0.11423, 0.2, 0.8];      % Specific mu values

% Preallocate for results
f_numeric_values = zeros(length(mua_range), length(mu_values));

% Loop through mu values and compute f_numeric
for i = 1:length(mu_values)
    mu_fixed = mu_values(i);  % Current mu value
    f_numeric_values(:, i) = f_numeric(mu_fixed, mua_range); % Compute for fixed mu
end

% Plot the results
figure;
hold on;
for i = 1:length(mu_values)
    plot(mua_range, f_numeric_values(:, i), 'DisplayName', ['\mu = ', num2str(mu_values(i))]);
end
hold off;

% Add labels, legend, and title
xlabel('\mu_a')
ylabel('fitness')
%title('2D Plots of \mu_a vs f_{numeric} at Different \mu Values')
legend show
grid on

%%
close all
% Define the range for mu and fixed values for mua
mu_range = linspace(0, 1, 500);    % Range for mu
mua_values = [0.00001, 0.004008, 0.02]; % Different fixed values of mua

% Initialize the figure
figure;
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
figure;
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
% close all
% % Define ranges for mu and mua
% mu_range = linspace(0, 1, 500);    % Range for mu
% mua_range = linspace(0,1, 500);  % Range for mua
% 
% % Create meshgrid
% [mu_mesh, mua_mesh] = meshgrid(mu_range, mua_range);
% 
% % Evaluate f_numeric on the meshgrid
% Z = f_numeric(mu_mesh, mua_mesh);
% 
% % Find the maximum value of Z and its location
% [max_Z, idx] = max(Z(:));
% [max_mu_idx, max_mua_idx] = ind2sub(size(Z), idx);
% max_mu = mu_mesh(max_mu_idx, max_mua_idx);
% max_mua = mua_mesh(max_mu_idx, max_mua_idx);
% 
% % Plot the contour plot
% figure;
% contour(mu_mesh, mua_mesh, Z, 20); % Contour plot with 20 levels
% hold on;
% 
% % Highlight the maximum value (argmax) with a red dot
% plot(max_mu, max_mua, 'r*', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% 
% % Add labels, title, and grid
% xlabel('\mu');
% ylabel('\mu_a');
% title('Contour Plot of Z with Maximum Highlighted');
% grid on;
% hold off;
% 
% % Display the maximum value and location
% disp(['Maximum Z occurs at \mu = ', num2str(max_mu), ' and \mu_a = ', num2str(max_mua)]);
% disp(['Maximum Z value is ', num2str(max_Z)]);
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
figure;
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