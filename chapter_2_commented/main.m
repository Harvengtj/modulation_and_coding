%% MODULATION AND CODING : MAIN (CHAPTER 2)
clc;clear;close all;

% --- Parameters ---

Nbps = 4;                                        % Number of bits per symbol
rollOff = 0.2;                                   % Roll-Off Factor
M = 100;                                           % Upsampling Factor
N = 16*M+1;                                            % Number of taps (ODD ONLY)
EbN0 = 10;                                 % Eb to N0 ratio 
symRate = 5e6;                               % Symbol Rate
Tsymb= 1/symRate;                        % Symbol Period
Fs = symRate*M;                                 % Sampling Frequency
Nsymb = 1000;                         % Number of symbols
Nb = Nbps*Nsymb;
Fc = 600e6;
ppm = 10;
CFO = ppm*Fc*1e-6;                              % Carrier Frequency Offset
CPE_deg = 0;
CPE_rad = CPE_deg*pi/180;
K = 0.05;
timeShift = 20;
Nb_iter = 50;
averagedTimeError = zeros(Nb_iter,Nsymb);
averagedBER = 0;

pilotPos = 34;
pilotLength = 20;
avgWindowLength = 8;

for iter = 1:Nb_iter

    %% TRANSMITTER
    % --- Bit Generation ---

    fprintf('Iteration n°%d\n',iter);
    bits_tx = randi(2,1,Nb)-1;               % bits_tx = Binary sequence

    
    % --- Mapping ---

    if Nbps > 1
            signal_tx = mapping(bits_tx.',Nbps,'qam').';         % Symbols sequence at transmitter
    else
            signal_tx = mapping(bits_tx.',Nbps,'pam').';         % Symbols sequence at transmitter   
    end
    

    % --- Divide message into unuseful data, pilot, symbols ---

    unuseful = signal_tx(1:pilotPos-1);
    pilot = signal_tx(pilotPos:pilotPos+pilotLength-1);
    symbols = signal_tx(pilotPos+pilotLength : end);


    % --- Upsampling ---

    upsampled_signal = zeros(1,length(signal_tx)*M);
    for i = 1:length(signal_tx)
        upsampled_signal(1+M*(i-1))=signal_tx(i);
        for j = 2:M
            upsampled_signal(j+M*(i-1))=0;
        end
    end

 
    % --- RRC Nyquist Filter TX ---

    [h_RRC,H_RRC] =  halfroot_Nyquist(Fs,Tsymb,N,rollOff);
    filtered_signal_tx = conv(upsampled_signal,h_RRC);

    %% AWGN channel
    % --- Noise ---

    avSymbEnergyBaseband = mean(abs(filtered_signal_tx).^2) * Tsymb;        % Calculate the average baseband symbol energy
    avSymbEnergy = (1/2) * avSymbEnergyBaseband;                            % Calculate the average symbol energy
    Eb = avSymbEnergy / Nbps;                                               % Calculate the energy per bit
    
    N0 = Eb ./ (10.^(EbN0/10));                                             % Calculate the noise PSD for each Eb/N0 value
    NoisePower = 2 * N0 * Fs;                                               % Calculate the noise power

    noise = sqrt(NoisePower/2).*(randn(1,length(signal_tx)*M+N-1)+1i*randn(1,length(signal_tx)*M+N-1));
    signal_rx = filtered_signal_tx + noise;

    % --- CFO and CPE ---

    t1 = ((0:size(signal_rx,2)-1))*1/Fs;
    signal_rx_sync_errors = signal_rx.*exp(1j*(2*pi*CFO.*t1+CPE_rad));

    %scatterplot(signal_rx_sync_errors);
    
    %% RECEIVER
    % --- RRC Nyquist Filter RX ---

    filtered_signal_rx = conv(signal_rx_sync_errors,fliplr(h_RRC));
    cropped_filtered_signal_rx = filtered_signal_rx(N:end-(N-1));
    

    % --- Time Shift ---
    
    downsampling_ratio=M/2;

    shifted_signal_rx=circshift(cropped_filtered_signal_rx,timeShift);
    partial_downsampled_signal_rx = downsample(shifted_signal_rx,downsampling_ratio);                                                  
    
    % --- Gardner ---

    [downsampled_signal_rx_corrected,time_error]=gardner(partial_downsampled_signal_rx,K,M/downsampling_ratio);
    averagedTimeError(iter,:) = time_error;
    

    % --- Data Acquisition ---

    [toa, est_CFO] = dataAcquisition(downsampled_signal_rx_corrected,pilot,avgWindowLength, Tsymb);
    disp(['ToA : ' num2str(toa)]);
    disp(['Estimated CFO : ' num2str(est_CFO)]);
    t2 = ((0:length(signal_tx)-1))*Tsymb;
    compensated_signal_rx = downsampled_signal_rx_corrected.*exp(-1j*2*pi*est_CFO*t2);
    scatterplot(compensated_signal_rx);
    

    % --- Demapping ---
    
    if Nbps > 1
        bits_rx = demapping(compensated_signal_rx.',Nbps,"qam");    
    else
        bits_rx = demapping(real(compensated_signal_rx.'),Nbps,"pam");
    end
    

    % --- BER ---

    BER = 0;
    
    for i = 1:Nb
        if(bits_rx(i) ~= bits_tx(i))
            BER = BER+1/Nb;
        end
    end
    
  
    averagedBER = averagedBER+BER/Nb_iter;
end
    disp(['BER = ' averagedBER]);


% --- Plot Error Convergence ---

MeanTimeError = mean(averagedTimeError);
MeanTimeError = (MeanTimeError + timeShift/M);

vector = 1:25:Nsymb;
plot(vector,MeanTimeError(vector),'ro-');
hold on;
legend(['CFO=' num2str(ppm) 'ppm']);

grid on;
legend('show');
xlabel("Symbols");
ylabel("Time error (mean)");
title('Time error VS Samples');
