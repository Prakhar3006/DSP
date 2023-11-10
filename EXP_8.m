%1
R0 = 2.5;
numerator = 1;
denominator = [1, -R0];
H1 = tf(numerator, denominator);
figure;
pzmap(H1);
title('Pole-Zero Plot');

%b
% Initialize variables
n = 20; % Number of days
y = zeros(1, n);
y(1) = 1; % Assuming initial infection on day 1

% Iterate to find daily infections
for k = 2:n
    y(k) = 1 + R0 * y(k - 1);
end

% Display the results
disp(y);

%c
threshold = 1e6; % 1 million
days_to_reach_threshold = ceil(log(threshold) / log(R0)) + 1;

disp(['It takes approximately ', num2str(days_to_reach_threshold), ' days to reach 1 million new daily infections.']);

%d
% one-point method
Y = 34285612;
R01 = 1-1/(Y)

%linear-regression


%e
n = 0:20;
y = R0.^n;

% Design an integrator filter
total_infections = cumsum(y);

% Plot the results
figure;
subplot(2, 1, 1);
stem(n, y, 'o-', 'LineWidth', 2);
xlabel('Days (n)');
ylabel('New Daily Infections');
title('Number of Newly Infected People Over Time');

subplot(2, 1, 2);
stem(n, total_infections, 'o-', 'LineWidth', 2);
xlabel('Days (n)');
ylabel('Total Infections');
title('Total Number of Infections Over Time');


% Problem 2 - Multi-pole IIR Filter

% Define the coefficients for the transfer function HM(z)
ak = [.1, .15, .25, .26, .34, .42, .25, .2, .15, .1, .1, .1];
M = length(ak);

% Part 1: Plot new daily infections for the first n = 100 days
n_days = 100;
n = 0:n_days;

% Kronecker delta as input
delta = (n == 0);

% Implement the IIR filter
y_multi_pole = filter(1, [1, -ak], delta);

% Part 2: Use an integrator filter to obtain total infections
filtered_output_multi_pole = filter(1, [1, -1], y_multi_pole);

% Plot new daily infections and total infections
figure;
subplot(2, 1, 1);
stem(n, y_multi_pole, 'b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Days');
ylabel('Number of New Daily Infections');
title('New Daily Infections (Multi-pole IIR Filter)');
grid on;

subplot(2, 1, 2);
plot(n, filtered_output_multi_pole, 'g--', 'LineWidth', 2);
xlabel('Days');
ylabel('Total Infections');
title('Total Infections (Integrated)');
grid on;





%3
% Define the coefficients and parameters
M = 12;
ak = [0.1, 0.15, 0.25, 0.26, 0.34, 0.42, 0.25, 0.2, 0.15, 0.1, 0.1, 0.1];
n_max = 100;
rho_values = [0.25, 0.50, 0.75];
total_infections = zeros(length(rho_values), 1);

% Initialize arrays for the input (Kronecker delta) and output (daily infections)
x = zeros(1, n_max + 1);
y = zeros(1, n_max + 1);

% Loop over different ρ values
for rho_idx = 1:length(rho_values)
    rho = rho_values(rho_idx);
    
    % Set the Kronecker delta at day 0
    x(1) = 1;
    
    % Apply the IIR filter with scaled coefficients
    for n = 1:n_max
        for k = 1:M
            if n - k > 0
                y(n + 1) = y(n + 1) + (1 - rho) * ak(k) * y(n - k + 1);
            end
        end
        y(n + 1) = 1 - y(n + 1);
    end
    
    % Calculate the total number of infections for n = 100 days
    total_infections(rho_idx) = sum(y);
    
    % Plot the daily infections for the current ρ value
    subplot(1, length(rho_values), rho_idx);
    plot(0:n_max, y);
    title(['ρ = ' num2str(rho)]);
    xlabel('Day');
    ylabel('Daily Infections');
    grid on;
end

% Display the total number of infections for each ρ value
disp('Total Infections for Different ρ Values:');
disp(total_infections);

%4
% Define parameters
R0 = 1.15;  % Reproduction number
K = 1e6;   % Population size
n_max = 100;  % Number of days

% Initialize arrays for the logistic and first-order models
x_logistic = zeros(1, n_max + 1);
x_first_order = zeros(1, n_max + 1);

% Initialize parameters for derivative calculation
D1_filter = [1, -1];  % First derivative filter coefficients
D2_filter = [1, -2, 1];  % Second derivative filter coefficients

% Apply the logistic model
for n = 0:n_max
    x_logistic(n + 1) = K / (1 + (K * (R0 - 1) - R0) * R0^(-(n + 1)) / (R0 - 1));
end

% Apply the first-order model
M = 12;  % Order of the first-order filter (adjust as needed)
ak = ones(1, M);  % Coefficients for the first-order model
x_first_order(1) = 1;  % Initial condition

for n = 1:n_max
    for k = 1:M
        if n - k > 0
            x_first_order(n + 1) = x_first_order(n + 1) + ak(k) * x_first_order(n - k + 1);
        end
    end
    x_first_order(n + 1) = 1 - x_first_order(n + 1);
end

% Plot the results
figure;
plot(0:n_max, x_logistic, 'b', 'LineWidth', 2);
hold on;
plot(0:n_max, x_first_order, 'r--', 'LineWidth', 2);
title('Total Infections: Logistic vs. First-Order Model');
legend('Logistic Model', 'First-Order Model');
xlabel('Day');
ylabel('Total Infections');
grid on;

% Calculate first derivative and second derivative
first_derivative = conv(x_logistic, D1_filter, 'valid');
second_derivative = conv(x_logistic, D2_filter, 'valid');

% Find the inflection point
[~, max_derivative_index] = max(first_derivative);
zero_crossing_index = find(diff(sign(second_derivative)) == 2, 1);

disp(['Inflection point (First Derivative Maximum): Day ' num2str(max_derivative_index)]);
disp(['Inflection point (Zero-Crossing of Second Derivative): Day ' num2str(zero_crossing_index)]);
