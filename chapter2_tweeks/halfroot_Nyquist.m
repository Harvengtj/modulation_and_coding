%% NYQUIST FILTER

function [h_RRC, H_RRC] = halfroot_Nyquist(Fs, Tsymb, N, RollOff)
% halfroot_Nyquist - Computes the Nyquist filter (RC) and its root (RRC)
%
% Inputs:
%   Fs      - Sampling frequency [Hz]
%   Tsymb   - Symbol duration [s]
%   N       - Number of FFT points (defines time resolution)
%   RollOff - Roll-off factor (0 < RollOff <= 1)
%
% Outputs:
%   h_RRC   - Time domain impulse response of the Raised Cosine filter
%   H_RRC   - Frequency domain representation of the Raised Cosine filter

%% (1) Frequency domain design
step_offset = (1/N)*Fs;
fmax = step_offset*((N-1)/2);
f = linspace(-fmax, fmax, N); % Frequency vector

% Initialize H_RC (Raised Cosine)
H_RC = zeros(1, N);

% Raised Cosine (RC) filter design
for i = 1:N
    abs_f = abs(f(i));
    if abs_f < (1 - RollOff) / (2 * Tsymb)
        H_RC(i) = Tsymb;
    elseif abs_f <= (1 + RollOff) / (2 * Tsymb)
        H_RC(i) = Tsymb * 0.5 * (1 + cos(pi * Tsymb * (abs_f - (1 - RollOff) / (2 * Tsymb)) / RollOff));
    else
        H_RC(i) = 0;
    end
end

%% (2) Time domain conversion for Raised Cosine filter
H_RC_plot = H_RC;
H_RC = ifftshift(H_RC);
H_RRC = sqrt(H_RC);       
h_RC = ifft(H_RC);
h_RRC = ifft(H_RRC);
h_RRC = fftshift(h_RRC/sqrt(h_RC(1)));
h_RC = fftshift(h_RC/h_RC(1)); 


%% (4) Plot impulse response and frequency response

% Time vector for plotting
Ts = 1 / Fs; % Sampling period

t = (-(N-1)/2:(N-1)/2) * Ts; % Time vector

t_symb = (-(N-1)/2:(N-1)/2) * Tsymb; % Time vector (s) if N even

% % Plot frequency response for the Raised Cosine filter
% figure;
% 
% 
% % Plot frequency response for H_RC (Raised Cosine) in magnitude (in dB)
% subplot(2,1,1);
% %plot(f, 20*log(abs(H_RC_plot)), 'b', 'LineWidth', 1.5); % Magnitude in dB (log scale)
% title('Frequency response of Raised Cosine filter');
% xlabel('Frequency [Hz]');
% ylabel('Magnitude [dB]');
% grid on;
% 
% % Plot the time domain impulse response for Raised Cosine filter
% subplot(2,1,2);
% %plot(t, h_RC, 'b', 'LineWidth', 1.5);
% hold on
% one_zeros = zeros(1, N);
% one_zeros((N+1)/2) = 1;
% stem(t_symb, one_zeros, 'r', MarkerSize=3, MarkerFaceColor='auto');
% hold on
% %plot(t, h_RRC, 'g', 'LineWidth', 1.5);
% title('Impulse response of Raised Cosine filter');
% 
% % Définir le texte de la textbox
% txt = {['N = ' num2str(N)], ['Roll-off = ' num2str(RollOff)]};
% 
% % Ajuster la position de la textbox (déplacer plus à gauche)
% annotation('textbox', [0.2, 0.7, 0.18, 0.15], 'String', txt, 'EdgeColor', 'none', 'BackgroundColor', 'white');
% 
% xlabel('Time [s]');
% ylabel('Normalized Amplitude');
% grid on;
% %xlim([-0.02 0.02]);
% 
% % Ajout des légendes
% legend('Raised Cosine (h_{RC})', 'Dirac Impulse', 'Root Raised Cosine (h_{RRC})');

end



