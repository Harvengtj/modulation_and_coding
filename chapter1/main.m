%% === MAIN ===
clc;clear;close all;


%% === Parameters ===
disp('=== Parameters ===');
Nbps = 4;  % Number of bits per symbol
order = 2^Nbps;  % Modulation order
Nsymb = 10000;  % Number of bits
Nb = Nbps*Nsymb;  % Number of bits
rollOff = 0.2;  % Roll-Off Factor
M = 4;  % Upsampling Factor
N = 101;  % Number of taps (ODD ONLY)
EbN0 = 0:2:16;  % Eb to N0 ratio (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
averageNb = 1;  % Number of iteration to average the BER   -> more than 1 to make an average
symbRate= 5e6;  % Symbol rate [symb/s]
Tsymb= 1/symbRate;  % Symbol Period
Fs = symbRate*M;  % Sampling Frequency


fprintf('Number of symbols : %d\nNumber of bits per symbols : %d [bit/symb]\nRoll-off factor : %d\nUpsampling factor : %d\nNumber of taps : %d\nSymbol rate : %d [symb/s]\n', Nsymb, Nbps, rollOff, M, N, symbRate);


% To display graphs, put averageNb=1, Nbps = int and EbN0 = int

% ============================================= %


%% === Bit Generation ===
% Generate random bits (Nbps bits per symbol)
bits_tx = randi([0 1], Nb, 1);
size(bits_tx)

%% === Mapping ===
fprintf('\n=== Mapping ===\n');
if Nbps > 1
    modulation = 'qam';
else
    modulation = 'pam';         
end

fprintf('Modulation type : %s\n', modulation);
signal_tx = mapping(bits_tx,Nbps,modulation); % Symbols sequence at transmitter
size(signal_tx)
scatterplot(signal_tx); % Constellation at transmitter side

%% === Upsampling === 
fprintf('\n=== Upsampling ===\n');
upsampled_signal_tx = upsample(signal_tx, M);

size(upsampled_signal_tx)
fprintf('Signal length : %d (= Nsymb*M)\n', length(upsampled_signal_tx));

%% === Nyquist Filter TX ===
fprintf('\n=== Nyquist Filter TX ===\n');
[h_RRC,H_RRC] =  halfroot_Nyquist(Fs, Tsymb, N, rollOff);
%[h_RRC,H_RRC] =  RRC(Fs, Tsymb, N, rollOff, Nbps, averageNb, M);
% h_RRC = sqrt(rcosdesign(rollOff, N, Fs/symbRate));
h_RRC = h_RRC';
filtered_signal_tx = conv(upsampled_signal_tx,h_RRC);
size(h_RRC)
size(filtered_signal_tx)
fprintf('Signal length : %d (= Nsymb*M + (N-1))\n', length(filtered_signal_tx));

%% === Additive White Gaussian Noise (AWGN) ===
fprintf('\n=== Additive White Gaussian Noise (AWGN) ===\n');
signalEnergyBaseband = (sum(abs(filtered_signal_tx).^2))*(1/Fs); % See "Representation of pass-band signals" : Es = sigma_s^2 * Ts (factor 1/2 because baseband --> RF)
signalEnergy = (1/2)*signalEnergyBaseband;
Eb = signalEnergy/(Nb);

N0 = Eb./(10.^(EbN0/10)); % <--> SNR = 10*log(Eb/N0)
noisePower = 2*N0*Fs;

noise = zeros(length(EbN0),Nsymb*M + (N-1)); % Because the halfroot Nyquist filter extends the signal with (N-1) extra elements
signal_rx = zeros(length(EbN0),Nsymb*M + (N-1));
size(signal_rx)

for k = 1:length(EbN0)
    noise(k,:) = sqrt(noisePower(k)/2).*(randn(1,Nsymb*M + (N-1)) + 1i*randn(1,Nsymb*M + (N-1)));
    signal_rx(k,:) = filtered_signal_tx' + noise(k,:);
end


%% === Nyquist Filter RX ===
fprintf('\n=== Nyquist Filter RX ===\n');
filtered_signal_rx = zeros(length(EbN0), Nsymb*M + 2*(N-1));
cropped_filtered_signal_rx = zeros(length(EbN0), Nsymb*M);
for i = 1:length(EbN0)
    filtered_signal_rx(i,:) = conv(signal_rx(i,:),fliplr(h_RRC'));
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(i,N:end-(N-1));
end

%% === Downsampling ===
fprintf('\n=== Downsampling ===\n');
downsampled_signal_rx = zeros(length(EbN0),Nsymb);
for j = 1:length(EbN0)
    for i = 1:Nsymb
        downsampled_signal_rx(j,i)=cropped_filtered_signal_rx(j,1 + M*(i-1));
    end
end

% Constellation at receiver side
for i = 1:length(EbN0)
    scatterplot(downsampled_signal_rx(i, :));
    title(['Constellation at Eb/N0 = ', num2str(EbN0(i)), ' dB']);
end 

%% === Demapping ===
downsampled_signal_rx = downsampled_signal_rx';
bits_rx = zeros(Nb, length(EbN0));
size(downsampled_signal_rx)
size(bits_rx)
for i = 1:length(EbN0)
    bits_rx(:, i) = demapping(downsampled_signal_rx(:, i), Nbps, modulation);
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

berEst = zeros(length(EbN0));

for i = 1:length(EbN0)
    % Calculate the number of bit errors
    nErrors = biterr(bits_tx,bits_rx(:, i));

    % Increment the error and bit counters
    numErrs = nErrors;
    numBits = Nb;

    % Estimate the BER
    berEst(i) = numErrs/numBits;

end


berTheory = berawgn(EbN0,modulation,order);

figure; % Ouvre une nouvelle figure
semilogy(EbN0, berEst, '*'); % Marqueurs rouges plus gros
hold on;
semilogy(EbN0, berTheory); % Ligne bleue plus épaisse
hold off;
grid on; % Active la grille

legend('Estimated BER', 'Theoretical BER');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
title('BER vs Eb/N0');




