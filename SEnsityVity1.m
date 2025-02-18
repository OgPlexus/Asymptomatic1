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
pri=linspace(1, 20, 100);
y_values = 1./pri;%linspace(0, 1, 100); % 100 points in [0, 1]
x_max_values = zeros(size(y_values)); % Preallocate for solutions
g_max_values = zeros(size(y_values)); % Preallocate for max g(x)

% Find argmax x for each y
for i = 1:length(y_values)
    % Substitute the current y value into g(x)
    g_y = subs(g, y, y_values(i));
    
    % Maximize g(x) over the interval [0, 1]
    [x_max, g_max] = fminbnd(@(x_eval) -double(subs(g_y, x, x_eval)), 0, 1);
    
    % Store the results
    x_max_values(i) = x_max;       % x that maximizes g(x)
    g_max_values(i) = -g_max;      % Maximum value of g(x)
end

% Plot y vs. x_max (argmax x)
%close all
figure;
plot(1./y_values, x_max_values, 'k-', 'LineWidth', 1.5);
xlabel('$\nu^{-1}$: mean period in symptomatic stage');
ylabel('virulence: ESS');
%title('x = argmax g(x) vs. y');
grid on;
hold on

% Optional: Plot y vs. max g(x)
% figure;
% plot(y_values, g_max_values, 'b-', 'LineWidth', 1.5);
% xlabel('y');
% ylabel('max g(x)');
% title('max g(x) vs. y');
% grid on;
% hold on