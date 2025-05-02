%% === MAIN ===
clc; clear; close all;

%% === Parameters ===
disp('=== Parameters ===');
Nbps = 4;               % Number of bits per symbol
order = 2^Nbps;         % Modulation order, calculated as 2^Nbps (e.g., 4 bits/symbol -> order 16)
Nsymb = 1000;          % Total number of symbols to transmit
Nb = Nbps * Nsymb;      % Total number of bits, calculated as Nbps * Nsymb
rollOff = 0.2;          % Roll-off factor for the Nyquist filter, determines the transition band
M = 100;                  % Upsampling factor, to increase the sampling frequency
N = 101;                % Number of taps in the Nyquist filter (must be odd)
bandwidth = 6e6;        % Cut-off frequency of the Nyquist filter
EbN0 = 20;          % Vector of Eb/N0 ratios (energy per bit to noise PSD) for BER comparison
symbRate = 5e6;         % Symbol rate [symb/s]
Tsymb = 1 / symbRate;   % Symbol period, inverse of the symbol rate
Fs = symbRate * M;      % Sampling frequency, calculated as symbRate * M
Fc = 2e9;             % Carrier frequency
ppm = [1*(1e-6) 10*(1e-6) 50*(1e-6) 100*(1e-6) 200*(1e-6)];         % PPM
phi_0 = 2*pi*rand; % Phase offset 
K = 0.05;               % Kappa
timeShift = 10;
nb_iter = 50;
AverageTimeError = zeros(nb_iter,length(ppm),Nsymb);


% Display parameters in the console for verification
fprintf('Number of symbols : %d\nNumber of bits per symbols : %d [bit/symb]\nRoll-off factor : %d\nUpsampling factor : %d\nNumber of taps : %d\nBandwidth : %d [Hz]\nSymbol rate : %d [symb/s]\n', Nsymb, Nbps, rollOff, M, N, bandwidth, symbRate);

% ============================================= %

for iter = 1:nb_iter

fprintf('\n---------- Iteration %d ----------\n',iter);
%% === Bit Generation ===
% Generate a random sequence of bits (Nb bits)
bits_tx = randi([0 1], Nb, 1);                                          % Use randi to generate bits 0 or 1

%% === Mapping ===
fprintf('\n=== Mapping ===\n');
if Nbps > 1
    modulation = 'qam';                                                 % Use QAM modulation if more than one bit per symbol
else
    modulation = 'pam';                                                 % Use PAM modulation if only one bit per symbol
end

fprintf('Modulation type : %s\n', modulation);                          % Display the modulation type
signal_tx = mapping_comments(bits_tx, Nbps, modulation);                % Map bits to symbols using the mapping function
%scatterplot(signal_tx);                                                 % Display the constellation diagram of the transmitted symbols


%% === Upsampling ===
fprintf('\n=== Upsampling ===\n');
upsampled_signal_tx = upsample(signal_tx, M);                           % Upsample the symbol sequence by a factor of M

%% === Nyquist Filter TX ===
fprintf('\n=== Nyquist Filter TX ===\n');
[h_RRC, H_RRC] = halfroot_Nyquist_comments(Fs, Tsymb, N, rollOff);      % Design a Root Raised Cosine (RRC) filter
h_RRC = h_RRC';                                                         % Transpose the filter coefficients for use in convolution
filtered_signal_tx = conv(upsampled_signal_tx, h_RRC);                  % Convolve the upsampled signal with the RRC filter

%% === Additive White Gaussian Noise (AWGN) ===
fprintf('\n=== Additive White Gaussian Noise (AWGN) ===\n');
avSymbEnergyBaseband = mean(abs(filtered_signal_tx).^2) * Tsymb;        % Calculate the average baseband symbol energy
avSymbEnergy = (1/2) * avSymbEnergyBaseband;                            % Calculate the average symbol energy
Eb = avSymbEnergy / Nbps;                                               % Calculate the energy per bit

N0 = Eb ./ (10.^(EbN0/10));                                             % Calculate the noise PSD for each Eb/N0 value
noisePower = 2 * N0 * Fs;                                               % Calculate the noise power

noise = sqrt(noisePower/2) .* (randn(1, Nsymb * M + (N-1)) + 1i * randn(1, Nsymb * M + (N-1))); % Generate AWGN                        % Initialize a matrix for the noise
signal_rx = filtered_signal_tx' + noise;  % Add noise to the filtered signal


%% === Carrier Frequency Offset (CFO) + Carrier Phase Error (CPE) ===
fprintf('\n=== Carrier Frequency Offset (CFO) + Carrier Phase Error (CPE) ===\n');

% Temporal axis
t_axis = (0:length(signal_rx) - 1)*(1/Fs); % Time axis


signal_rx_sync_errors = zeros(length(ppm), length(signal_rx));

for k=1:length(ppm)
    % CFO
    CFO = exp(1i*2*pi*(ppm(k)*Fc)*t_axis); 
    % CPE 
    CPE = exp(1i*phi_0); 
    
    % Total shift
    shift = CFO .* CPE;
    signal_rx_sync_errors(k,:) = signal_rx .* shift;
end



%% === Nyquist Filter RX ===
fprintf('\n=== Nyquist Filter RX ===\n');
filtered_signal_rx = zeros(1, Nsymb * M + 2 * (N-1));        % Initialize a matrix for the filtered signal
cropped_filtered_signal_rx = zeros(length(ppm), Nsymb * M);            % Initialize a matrix for the cropped signal
for i = 1:length(ppm)
    filtered_signal_rx = conv(signal_rx_sync_errors(i,:), fliplr(h_RRC'));     % Convolve the received signal with the flipped RRC filter
    cropped_filtered_signal_rx(i,:) = filtered_signal_rx(N:end-(N-1)); % Crop the filtered signal to remove extensions (causality)
end

    %% === Time shift ===
    downsampling_ratio = M/2;
    shifted_signal_rx = zeros(length(ppm),length(signal_tx)*M);
    partial_downsampled_signal_rx = zeros(length(ppm),length(signal_tx)*M/downsampling_ratio);  % keep 2 samples per symbol (required for Gardner)

    for k = 1:length(ppm)
            shifted_signal_rx(k,:) = circshift(cropped_filtered_signal_rx(k,:),timeShift);
            partial_downsampled_signal_rx(k,:) = downsample(shifted_signal_rx(k,:),downsampling_ratio); % downsample by M/2
    end                                                          
    
    
   
    %% === Gardner ===
    downsampled_signal_rx_corrected = zeros(length(ppm),length(signal_tx));
    time_error = zeros(length(ppm),length(signal_tx));
    for k = 1:length(ppm)
        [downsampled_signal_rx_corrected(k,:),time_error(k,:)] = gardner(partial_downsampled_signal_rx(k,:),K,M/downsampling_ratio);
        AverageTimeError(iter,k,:) = time_error(k,:);
    end
    
end


% Compute mean and standard deviation
MeanTimeError = zeros(length(ppm),Nsymb);
VarianceTimeError = zeros(length(ppm),Nsymb);
    for k = 1:length(ppm)
        MeanTimeError(k,:) = mean(AverageTimeError(:,k,:));
        MeanTimeError(k,:) = (MeanTimeError(k,:)+timeShift/M);
        VarianceTimeError(k,:) = std(AverageTimeError(:,k,:));
    end
colorVector = ['r','b','g','m','c','k','y'];
    Legend=cell(length(ppm));

% Plot mean and standard deviation curves
for i = 1:length(ppm)
   vector = 1:25:Nsymb;

   p(i) = plot(vector,MeanTimeError(i, vector),[colorVector(i) 'o-']);
   hold on;
   plot(vector,MeanTimeError(i,vector)-VarianceTimeError(i,vector),[colorVector(i) '--']);
   hold on;
   plot(vector,MeanTimeError(i,vector)+VarianceTimeError(i,vector),[colorVector(i) '--']);
   hold on;
   Legend{i}=['CFO=' num2str(ppm(i))];
end

legend(p(1:length(ppm)),Legend(1:length(ppm)));

grid on;

legend('show');
xlabel("Symbols");
ylabel("Time error (mean±deviation)");
if(Nbps==1) 
    text='BPSK ';
elseif(Nbps==2) 
    text='QPSK ';
elseif(Nbps==4) 
    text='16QAM ';
else
    text ='64QAM ';
end
% 
% txt = {['#taps= ' num2str(N)],['RollOff= ' num2str(RollOff)],['M= ' num2str(M)],['SymRate= ' num2str(SymRate*1e-6) 'MBd'],['Phase Offset=' num2str(phase_offset_deg) '°']};
% annotation('textbox',[0.2,0.2,0.22,0.22],'String',txt,'BackgroundColor','white');

title([text,'(Nbps=',num2str(Nbps),')']);