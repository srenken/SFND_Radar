% Parameters
N_train = 10;  % Number of training cells on each side of the CUT
N_guard = 2;   % Number of guard cells on each side of the CUT
P_fa = 1e-3;   % Probability of false alarm
Ns = 1000;     % Number of samples in the signal

% Generate random noise + target signal
noise = abs(randn(1, Ns));
signal = noise;  % Initialize signal as noise
target_position = 500;
signal(target_position) = signal(target_position) + 8;  % Insert target

% Number of cells in the window (training + guard cells on both sides)
N_cells_total = 2 * (N_train + N_guard) + 1;

% Threshold factor based on the desired P_fa
threshold_factor = (N_train * 2) * (P_fa^(-1 / (2 * N_train)) - 1);

% Initialize detection results
CFAR_output = zeros(1, Ns);  % Store detection results

% Slide window across the signal from index 1 to (Ns - N_train - N_guard)
for CUT = 1:(Ns - (N_guard + N_train + 1))
    
    % Check if the CUT is within valid range for training and guard cells
    if CUT <= N_train + N_guard || CUT > Ns - (N_train + N_guard)
        continue;  % Skip invalid CUT positions
    end
    
    % Define the training cells
    leading_train_cells = signal(max(1, CUT - N_train - N_guard) : CUT - N_guard - 1);
    lagging_train_cells = signal(CUT + N_guard + 1 : min(Ns, CUT + N_guard + N_train));
    
    % Average noise level in the training cells
    noise_level = mean([leading_train_cells, lagging_train_cells]);
    
    % Calculate the threshold
    threshold = threshold_factor * noise_level;
    
    % Compare CUT with the threshold
    if signal(CUT) > threshold
        CFAR_output(CUT) = 1;  % Detection occurs
    end
end

% Plot the results
figure;
subplot(2,1,1);
plot(signal);
title('Signal with Target');
xlabel('Cell Index');
ylabel('Amplitude');

subplot(2,1,2);
plot(CFAR_output);
title('CFAR Detection Results');
xlabel('Cell Index');
ylabel('Detection (1 = Detected, 0 = Not Detected)');