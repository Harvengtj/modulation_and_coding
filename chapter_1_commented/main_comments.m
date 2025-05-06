%% === MAIN ===
clc; clear; close all;

%% === Parameters ===
disp('=== Parameters ===');
Nbps = 4;               % Number of bits per symbol
order = 2^Nbps;         % Modulation order, calculated as 2^Nbps (e.g., 4 bits/symbol -> order 16)
Nsymb = 10000;          % Total number of symbols to transmit
Nb = Nbps * Nsymb;      % Total number of bits, calculated as Nbps * Nsymb
rollOff = 0.2;          % Roll-off factor for the Nyquist filter, determines the transition band
M = 4;                  % Upsampling factor, to increase the sampling frequency
N = 101;                % Number of taps in the Nyquist filter (must be odd)
bandwidth = 6e6;        % Cut-off frequency of the Nyquist filter
EbN0 = 0:1:20;          % Vector of Eb/N0 ratios (energy per bit to noise PSD) for BER comparison
averageNb = 1;          % Number of iterations to average the BER
symbRate = 5e6;         % Symbol rate [symb/s]
Tsymb = 1 / symbRate;   % Symbol period, inverse of the symbol rate
Fs = symbRate * M;      % Sampling frequency, calculated as symbRate * M

% Display parameters in the console for verification
fprintf('Number of symbols : %d\nNumber of bits per symbols : %d [bit/symb]\nRoll-off factor : %d\nUpsampling factor : %d\nNumber of taps : %d\nBandwidth : %d [Hz]\nSymbol rate : %d [symb/s]\n', Nsymb, Nbps, rollOff, M, N, bandwidth, symbRate);

% ============================================= %

%% === Bit Generation ===
% Generate a random sequence of bits (Nb bits)
bits_tx = randi([0 1], Nb, 1);                                          % Use randi to generate bits 0 or 1
size(bits_tx)                                                           % Display the size of the generated bit sequence

%% === Mapping ===
fprintf('\n=== Mapping ===\n');
if Nbps > 1
    modulation = 'qam';                                                 % Use QAM modulation if more than one bit per symbol
else
    modulation = 'pam';                                                 % Use PAM modulation if only one bit per symbol
end

fprintf('Modulation type : %s\n', modulation);                          % Display the modulation type
signal_tx = mapping_comments(bits_tx, Nbps, modulation);                         % Map bits to symbols using the mapping function
size(signal_tx)                                                         % Display the size of the symbol sequence
scatterplot(signal_tx);                                                 % Display the constellation diagram of the transmitted symbols

%% === Upsampling ===
fprintf('\n=== Upsampling ===\n');
upsampled_signal_tx = upsample(signal_tx, M);                           % Upsample the symbol sequence by a factor of M

size(upsampled_signal_tx)                                               % Display the size of the upsampled signal
fprintf('Signal length : %d (= Nsymb*M)\n', length(upsampled_signal_tx)); % Length of the signal after upsampling

%% === Nyquist Filter TX ===
fprintf('\n=== Nyquist Filter TX ===\n');
[h_RRC, H_RRC] = halfroot_Nyquist_comments(Fs, Tsymb, N, rollOff);               % Design a Root Raised Cosine (RRC) filter
h_RRC = h_RRC';                                                         % Transpose the filter coefficients for use in convolution
filtered_signal_tx = conv(upsampled_signal_tx, h_RRC);                  % Convolve the upsampled signal with the RRC filter
size(h_RRC)                                                             % Display the size of the filter coefficients
size(filtered_signal_tx)                                                % Display the size of the filtered signal
fprintf('Signal length : %d (= Nsymb*M + (N-1))\n', length(filtered_signal_tx)); % Length of the signal after filtering

%% === Additive White Gaussian Noise (AWGN) ===
fprintf('\n=== Additive White Gaussian Noise (AWGN) ===\n');
avSymbEnergyBaseband = mean(abs(filtered_signal_tx).^2) * Tsymb;        % Calculate the average baseband symbol energy
avSymbEnergy = (1/2) * avSymbEnergyBaseband;                            % Calculate the average symbol energy
Eb = avSymbEnergy / Nbps;                                               % Calculate the energy per bit

N0 = Eb ./ (10.^(EbN0/10));                                             % Calculate the noise PSD for each Eb/N0 value
noisePower = 2 * N0 * Fs;                                               % Calculate the noise power

noise = zeros(length(EbN0), Nsymb * M + (N-1));                         % Initialize a matrix for the noise
signal_rx = zeros(length(EbN0), Nsymb * M + (N-1));                     % Initialize a matrix for the received signal
size(signal_rx)                                                         % Display the size of the received signal matrix

for k = 1:length(EbN0)
    noise(k,:) = sqrt(noisePower(k)/2) .* (randn(1, Nsymb * M + (N-1)) + 1i * randn(1, Nsymb * M + (N-1))); % Generate AWGN
    signal_rx(k,:) = filtered_signal_tx' + noise(k,:);                  % Add noise to the filtered signal
end

%% === Nyquist Filter RX ===
fprintf('\n=== Nyquist Filter RX ===\n');
filtered_signal_rx = zeros(length(EbN0), Nsymb * M + 2 * (N-1));        % Initialize a matrix for the filtered signal
cropped_filtered_signal_rx = zeros(length(EbN0), Nsymb * M);            % Initialize a matrix for the cropped signal
for i = 1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:), fliplr(h_RRC'));     % Convolve the received signal with the flipped RRC filter
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i, N:end-(N-1)); % Crop the filtered signal to remove extensions
end

%% === Downsampling ===
fprintf('\n=== Downsampling ===\n');
downsampled_signal_rx = zeros(length(EbN0), Nsymb);                     % Initialize a matrix for the downsampled signal
for j = 1:length(EbN0)
    for i = 1:Nsymb
        downsampled_signal_rx(j, i) = cropped_filtered_signal_rx(j, 1 + M * (i-1)); % Downsample the signal
    end
end

% Constellation at receiver side
for i = 1:length(EbN0)
    scatterplot(downsampled_signal_rx(i, :));
    title(['Constellation at Eb/N0 = ', num2str(EbN0(i)), ' dB']);
end 

%% === Demapping ===
downsampled_signal_rx = downsampled_signal_rx';                         % Transpose the downsampled signal matrix
bits_rx = zeros(Nb, length(EbN0));                                      % Initialize a matrix for the received bits
size(downsampled_signal_rx)                                             % Display the size of the downsampled signal matrix
size(bits_rx)                                                           % Display the size of the received bits matrix
for i = 1:length(EbN0)
    bits_rx(:, i) = demapping_comments(downsampled_signal_rx(:, i), Nbps, modulation); % Demap symbols to bits
end

%% === Verification ===
for i = 1:length(EbN0)
    if isequal(bits_tx, bits_rx(:, i))
        disp('SUCCESS: The transmitted and received bit streams match.');
    else
        disp('ERROR: The transmitted and received bit streams do not match.');
    end
end

%% === BER Computation ===
fprintf('=== BER Computation ===\n');

berEst = zeros(length(EbN0));                                           % Initialize a vector for the estimated BER

for i = 1:length(EbN0)
    nErrors = biterr(bits_tx,bits_rx(:, i));                            % Calculate the number of bit errors
    numErrs = nErrors;                                                  % Increment the error and bit counters
    numBits = Nb;

    berEst(i) = numErrs/numBits;                                        % Estimate the BER

end

berTheory = berawgn(EbN0, modulation, order);                           % Calculate the theoretical BER

figure;
semilogy(EbN0, berEst, '*');                                            % Plot the estimated BER
hold on;
semilogy(EbN0, berTheory);                                              % Plot the theoretical BER
hold off;
grid on;
ylim([10^(-6) 1]);                                                      % Limit the y-axis for better visualization
legend('Estimated BER', 'Theoretical BER');                             % Add a legend
xlabel('Eb/N0 (dB)');                                                   % Label the x-axis
ylabel('Bit Error Rate');                                               % Label the y-axis
title('BER vs Eb/N0');                                                  % Add a title to the plot
