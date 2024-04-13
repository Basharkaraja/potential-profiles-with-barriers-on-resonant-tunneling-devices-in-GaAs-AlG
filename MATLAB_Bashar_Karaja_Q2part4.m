mystartdefaults
MATLAB_Bashar_Karaja_Q2part1

% Define the x array
x = linspace(x_min, x_max, n);
% Bias values
U1_values = [-0.2, 0.2];

% Initialize arrays to store the modified potentials
W_minus_02 = zeros(1, n);
W_plus_02 = zeros(1, n);

% Loop over the bias values
for index = 1:length(U1_values)
    U1 = U1_values(index);

    % Potential barrier
    U = zeros(1, n);
    U(x > 0 & x <= 15) = 0.2;
    U(x >= 65 & x <= 80) = 0.2;

    % Calculate the electric field and the modified potential
    E = -U1 / (x_prime_max - x_prime_min);
    W = U;
    W(x > x_prime_min & x < x_prime_max) = U(x > x_prime_min & x < x_prime_max) - E * x(x > x_prime_min & x < x_prime_max);

    % Store the modified potential for each bias
    if U1 == -0.2
        W_minus_02 = W;
    elseif U1 == 0.2
        W_plus_02 = W;
    end
end

for i = 1:length(W_minus_02)
    if x(i) >= 80
        W_minus_02(i) = -0.2;
    end
end

% Plot the modified potential profile for U1 = -0.2 eV
figure;
plot(x, W_minus_02, 'LineWidth', 3);
title('U1 = -0.2 eV');
xlabel('x(Å)');
ylabel('E(eV)');
ylim([-0.25 0.25]);

for i = 1:length(W_plus_02)
    if x(i) >= 80
        W_plus_02(i) = 0.2;
    end
end
% Plot the modified potential profile for U1 = 0.2 eV
figure;
plot(x, W_plus_02, 'LineWidth', 3);
title('U1 = 0.2 eV');
ylabel('E(eV)');
xlabel('x(Å)');
ylim([-0.05 0.45]);