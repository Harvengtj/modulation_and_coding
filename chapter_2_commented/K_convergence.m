%% MODULATION AND CODING : KAPPA CONVERGENCE
clc;clear;close all;

%------Parameters------%
Nbps=4;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
RollOff= 0.2;                                   % Roll-Off Factor
M= 100;                                          % Upsampling Factor
N = 16*M + 1;                                     % Number of taps (ODD ONLY)
EbN0 = 1000;                                    % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
SymRate= 5e6;                               % Symbol Rate
Tsymb= 1/SymRate;                        % Symbol Period
Fs = SymRate*M;                                 % Sampling Frequency
Nsymb = 1000;
Nb = Nsymb*Nbps;
timeShift = 20;
Fc = 600e6;
ppm = 0;
CFO = ppm*Fc*1e-6;
phase_offset_deg = 0;
phase_offset= phase_offset_deg*pi/180;
K=[0.01 0.05];
AverageNb= 100;
AverageTimeError = zeros(AverageNb,Nb/Nbps,length(K));


for avr = 1:AverageNb
    %%
    % Bit Generation
    %------------------------
    
    fprintf('Iteration n°%d\n',avr);
    bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

    %%
    % Mapping
    %------------------------

    if Nbps>1
            signal_tx = mapping(bits_tx.',Nbps,'qam').';         % Symbols sequence at transmitter
    else
            signal_tx = mapping(bits_tx.',Nbps,'pam').';         % Symbols sequence at transmitter   
    end

    %%
    % Upsampling
    %-----------------

    upsampled_signal = zeros(1,length(signal_tx)*M);
    for i = 1:length(signal_tx)
        upsampled_signal(1+M*(i-1))=signal_tx(i);
        for j = 2:M
            upsampled_signal(j+M*(i-1))=0;
        end
    end

    %%
    % RRC Nyquist Filter TX
    %-------------------------

    [h_RRC,H_RRC] =  halfroot_Nyquist(Fs,Tsymb,N,RollOff);
    filtered_signal_tx = conv(upsampled_signal,h_RRC);

    %%
    % Noise
    %-----------------

    avSymbEnergyBaseband = mean(abs(filtered_signal_tx).^2) * Tsymb;        % Calculate the average baseband symbol energy
    avSymbEnergy = (1/2) * avSymbEnergyBaseband;                            % Calculate the average symbol energy
    Eb = avSymbEnergy / Nbps;                                               % Calculate the energy per bit
    
    N0 = Eb ./ (10.^(EbN0/10));                                             % Calculate the noise PSD for each Eb/N0 value
    NoisePower = 2 * N0 * Fs;                                               % Calculate the noise power

    noise = sqrt(NoisePower/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
    signal_rx = filtered_signal_tx + noise;
    
    %%
    % CFO & Carrier Phase Error
    %--------------------

    t1 = ((0:size(signal_rx,2)-1))*1/Fs;
    signal_rx_sync_errors = signal_rx.*exp(1j*(2*pi*CFO.*t1+phase_offset));

    %%
    % RRC Nyquist Filter RX
    %-------------------------

    t2=((0:length(signal_tx)*M-1))*1/Fs;
    filtered_signal_rx = conv(signal_rx_sync_errors,h_RRC);

    %%
    % Time Shift
    %-----------------------
    cropped_filtered_signal_rx = filtered_signal_rx(N:end-(N-1));
    shifted_signal_rx=circshift(cropped_filtered_signal_rx,timeShift);
    
    %cropped_filtered_signal_rx = cropped_filtered_signal_rx.*exp(-1j*(2*pi*CFO.*t2));
                                                                  
    
    
    %%
    % Gardner
    %-------------
    downsampling_ratio=M/2;
    partial_downsampled_signal_rx = downsample(shifted_signal_rx,downsampling_ratio);
    downsampled_signal_rx_corrected = zeros(length(signal_tx),length(K));
    time_error = zeros(length(signal_tx),length(K));
    for k = 1:length(K)
        [downsampled_signal_rx_corrected(:,k),time_error(:,k)]=gardner(partial_downsampled_signal_rx,K(k),M/downsampling_ratio);
        AverageTimeError(avr,:,k)=time_error(:,k);
    end

end
    MeanTimeError = zeros(Nb/Nbps,length(K));
    VarianceTimeError = zeros(Nb/Nbps,length(K));
    for k = 1:length(K)
        MeanTimeError(:,k) = mean(AverageTimeError(:,:,k));
        MeanTimeError(:,k)=(MeanTimeError(:,k) + timeShift/M);
        VarianceTimeError(:,k) = std(AverageTimeError(:,:,k));
    end
colorVector = ['k','r','b','g','y','m','c'];
    Legend=cell(length(K));
for i = 1:length(K)
   vector = 1:25:Nb/Nbps;

   p(i) = plot(vector,MeanTimeError(vector,i),[colorVector(i) '*-']);
   hold on;
   plot(vector,MeanTimeError(vector,i)-VarianceTimeError(vector,i),[colorVector(i) '--']);
   hold on;
   plot(vector,MeanTimeError(vector,i)+VarianceTimeError(vector,i),[colorVector(i) '--']);
   Legend{i}=['\kappa = ' num2str(K(i))];
end

legend(p(1:length(K)),Legend(1:length(K)));

grid on;
xlabel("Symbols");
ylabel("Time error (\mu ± \sigma)");
title('Time error std VS Sample index for different \kappa');