%% MODULATION AND CODING : ToA error VS SNR (including CFO)
clc;clear;close all;

%------Parameters------%
Nbps = 4;                                        % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6) -> vector to compare 
RollOff= 0.2;                                   % Roll-Off Factor
M = 100;                                          % Upsampling Factor
N = 16*M+1;                                     % Number of taps (ODD ONLY)
EbN0 = -6:2:20;                                    % Eb to N0 ratio  (Eb = bit energy, N0 = noise PSD)  -> vector to compare BER
SymRate= 5e6 ;                              % Symbol Rate
Tsymb= 1/(SymRate);                        % Symbol Period
Fs = SymRate*M;                                 % Sampling Frequency
Nsymb = 1000;          % Total number of symbols to transmit
Nb = Nbps * Nsymb;      % Total number of bits, calculated as Nbps * Nsymb
timeShift = [0 2];
Fc = 600e6;
ppm = [0 10];
CFO = ppm*Fc*1e-6;
phase_offset= 2*pi*rand;
K=0.05;
AverageNb= 50;


ToA = 34;
pilot_size = 20;
avgWindow_size = 8;
AverageTimeError = zeros(AverageNb,length(EbN0),length(CFO));

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
    %Divide msg into unuseful data, pilot, symbols
    %----------------------------------------------

    unuseful = signal_tx(1:ToA-1);
   
    pilot = signal_tx(ToA:ToA+pilot_size-1);
    symbols = signal_tx(ToA+pilot_size:end); % Actually not used
    
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

    noise = zeros(length(EbN0),length(signal_tx)*M+N-1);
    signal_rx = zeros(length(EbN0),length(signal_tx)*M+N-1);
    
    for i = 1:length(EbN0)
        noise(i,:) = sqrt(NoisePower(i)/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
        signal_rx(i,:) = filtered_signal_tx + noise(i,:);
    end
    
    %%
    % CFO & Carrier Phase Error
    %--------------------

    t1 = ((0:size(signal_rx,2)-1))*1/Fs;
    signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(CFO));
    for i=1:length(EbN0)
        for k=1:length(CFO)
           signal_rx_sync_errors(i,:,k) = signal_rx(i,:).*exp(1j*(2*pi*CFO(k).*t1+phase_offset));
        end
    end
     
    %%
    % RRC Nyquist Filter RX
    %-------------------------

    filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M+2*(N-1),length(CFO));
    for i =1:length(EbN0)
        for k =1:length(CFO)
            filtered_signal_rx(i,:,k) = conv(signal_rx_sync_errors(i,:,k),fliplr(h_RRC));
        end
    end                                                                      

    %%
    % Time Shift
    %-----------------------
    shifted_signal_rx = zeros(length(EbN0),length(filtered_signal_rx),length(CFO));
    cropped_filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*M,length(CFO));
    t2=((0:length(signal_tx)*M-1))*1/Fs;
    for i = 1:length(EbN0)
        for k =1:length(CFO)
            shifted_signal_rx(i,:,k)=circshift(filtered_signal_rx(i,:,k),timeShift(k));
            cropped_filtered_signal_rx(i,:,k) = shifted_signal_rx(i,N:end-(N-1),k);
            %cropped_filtered_signal_rx(i,:,k) = cropped_filtered_signal_rx(i,:,k).*exp(-1j*(2*pi*CFO.*t2));
        end
    end

%     %%
%     % Gardner
%     %-------------
%     downsampling_ratio=M/2;
%     partial_downsampled_signal_rx = zeros(length(EbN0),size(cropped_filtered_signal_rx,2)/downsampling_ratio);
%     downsampled_signal_rx_corrected = zeros(length(EbN0),length(signal_tx));
%     time_error = zeros(length(EbN0),length(signal_tx));
%     for i = 1:length(EbN0)
%         partial_downsampled_signal_rx(i,:) = downsample(cropped_filtered_signal_rx(i,:),downsampling_ratio);
%         [downsampled_signal_rx_corrected(i,:),time_error(i,:)]=gardner(partial_downsampled_signal_rx(i,:),K,M/downsampling_ratio);
%     end
    
    %%
    % Downsample
    %-----------------
    
    downsampled_signal_rx = zeros(length(EbN0),length(signal_tx),length(CFO));
    for i = 1:length(EbN0)
        for k =1:length(CFO)
            downsampled_signal_rx(i,:,k) = downsample(cropped_filtered_signal_rx(i,:,k),M);
        end
    end
    
    %%
    % Data acquisition
    %-----------------------------------

    for i = 1:length(EbN0)
        for k = 1:length(CFO)
            [est_ToA, est_CFO] = dataAcquisition(downsampled_signal_rx(i,:,k),pilot,avgWindow_size, Tsymb);
            AverageTimeError(avr,i,k) = abs(est_ToA-ToA);
        end
    end

    
end

VarianceTimeError = zeros(length(EbN0),length(CFO));
for i=1:length(CFO)
    for j = 1:length(EbN0)
        matrix = AverageTimeError(:,j,i).';
        VarianceTimeError(j,i) = std(matrix);
    end
end 

plot(EbN0,VarianceTimeError(:,1),'kx-');
hold on;
plot(EbN0,VarianceTimeError(:,2),'r<-');
hold off;
grid on;
xlabel("EbN0 [dB]");
ylabel("Time error stdev [samples]");%±deviation
legend('no CFO, \epsilon = 0', 'CFO = 10 ppm, \epsilon = 0.02', 'Location','northeast');
title('ToA error standard deviation VS SNR');