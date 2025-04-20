function [bit_rx] = demapping(symb_rx, Nbps, modulation)

% INPUTS:
% - symb_rx : vector of input symbols (variance 1)
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - bit_rx : vector of output bits

Nsymb = size(symb_rx, 1); % Calculate the number of symbols

switch modulation

    case 'pam'

        % Symbol to Integer Conversion
        sigma = sqrt(sum(([0:2^Nbps-1] - (2^Nbps-1)/2).^2) / 2^Nbps);   % Calculate the standard deviation
        int_rx = sigma * symb_rx + (2^Nbps-1)/2;                        % Scale and shift the symbols to convert them to integer values

        % Integer Detection
        int_det = round(int_rx);                                        % Round the scaled symbols to the nearest integer
        int_det(find(int_det < 0)) = 0;                                 % Clip any negative values to 0
        int_det(find(int_det > 2^Nbps-1)) = 2^Nbps-1;                   % Clip any values exceeding the maximum integer value

        % Integer to Binary Conversion
        mapp_rx = fliplr(de2bi(int_det));                               % Convert the integer values to binary, fliplr to reverse column order

        % Binary to Gray Conversion
        bit_rx2(:, 1) = mapp_rx(:, 1);                                  % The first column remains the same
        for ii = 2:Nbps
            bit_rx2(:, ii) = xor(mapp_rx(:, ii-1), mapp_rx(:, ii));     % XOR operation to convert binary to Gray code
        end

        bit_rx = reshape(bit_rx2', Nsymb * Nbps, 1);                    % Reshape the matrix into a vector of bits

    case 'qam'

        % REAL PART
        NbpsI = Nbps / 2;                                               % Number of bits per symbol for the real part
        symb_rxI = real(symb_rx);                                       % Extract the real part of the complex symbols

        % Symbol to Integer Conversion for the Real Part
        sigmaI = sqrt(sum(([0:2^NbpsI-1] - (2^NbpsI-1)/2).^2) / 2^NbpsI); % Calculate the standard deviation
        int_rxI = sigmaI * sqrt(2) * symb_rxI + (2^NbpsI-1)/2;          % Scale and shift the real part

        % Integer Detection for the Real Part
        int_detI = round(int_rxI);                                      % Round to the nearest integer
        int_detI(find(int_detI < 0)) = 0;                               % Clip negative values to 0
        int_detI(find(int_detI > 2^NbpsI-1)) = 2^NbpsI-1;               % Clip values exceeding the maximum integer

        % Integer to Binary Conversion for the Real Part
        mapp_rxI = fliplr(de2bi(int_detI));                             % Convert the integer values to binary, fliplr to reverse column order

        % Binary to Gray Conversion for the Real Part
        bit_rx2I(:, 1) = mapp_rxI(:, 1);                                % The first column remains the same
        for ii = 2:NbpsI
            bit_rx2I(:, ii) = xor(mapp_rxI(:, ii-1), mapp_rxI(:, ii));  % XOR operation to convert binary to Gray code
        end

        % IMAGINARY PART
        NbpsQ = Nbps / 2;                                               % Number of bits per symbol for the imaginary part
        symb_rxQ = imag(symb_rx);                                       % Extract the imaginary part of the complex symbols

        % Symbol to Integer Conversion for the Imaginary Part
        sigmaQ = sqrt(sum(([0:2^NbpsQ-1] - (2^NbpsQ-1)/2).^2) / 2^NbpsQ); % Calculate the standard deviation
        int_rxQ = sigmaQ * sqrt(2) * symb_rxQ + (2^NbpsQ-1)/2;          % Scale and shift the imaginary part

        % Integer Detection for the Imaginary Part
        int_detQ = round(int_rxQ);                                      % Round to the nearest integer
        int_detQ(find(int_detQ < 0)) = 0;                               % Clip negative values to 0
        int_detQ(find(int_detQ > 2^NbpsQ-1)) = 2^NbpsQ-1;               % Clip values exceeding the maximum integer

        % Integer to Binary Conversion for the Imaginary Part
        mapp_rxQ = fliplr(de2bi(int_detQ));                             % Convert the integer values to binary, fliplr to reverse column order

        % Binary to Gray Conversion for the Imaginary Part
        bit_rx2Q(:, 1) = mapp_rxQ(:, 1);                                % The first column remains the same
        for ii = 2:NbpsQ
            bit_rx2Q(:, ii) = xor(mapp_rxQ(:, ii-1), mapp_rxQ(:, ii));  % XOR operation to convert binary to Gray code
        end

        % BIT CONCATENATION
        bit_rx = reshape([bit_rx2I, bit_rx2Q]', Nsymb * Nbps, 1);       % Concatenate real and imaginary bits and reshape into a vector

end
