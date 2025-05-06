function [symb_tx] = mapping(bit_tx, Nbps, modulation)

% INPUTS:
% - bit_tx : vector of input bits
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - symb_tx : vector of output symbols (variance 1)

%% STEP 1
Nsymb = size(bit_tx, 1) / Nbps;             % Calculate the number of symbols

bit_tx2 = reshape(bit_tx, Nbps, Nsymb)';    % Reshape the bit vector into a matrix with Nbps rows and Nsymb columns, then transpose

%% EXAMPLE 1
% If bit_tx = [1; 0; 1; 1; 0; 0; 1; 0] and Nbps = 2, then:
%
% Nsymb = 8 / 2 = 4
% reshape(bit_tx, 2, 4) gives:
%
% 1  1  0  1
% 0  1  0  0
%
% Transposing it results in:
%
% 1  0
% 1  1
% 0  0
% 1  0

%% STEP 2
switch modulation

    case 'pam'

        % Gray to Binary Conversion
        mapp_tx(:, 1) = bit_tx2(:, 1);                                  % The first column remains the same
        for ii = 2:Nbps
            mapp_tx(:, ii) = xor(mapp_tx(:, ii-1), bit_tx2(:, ii));     % XOR operation to convert Gray code to binary
        end

        % Binary to Integer Conversion
        int_tx = bi2de(fliplr(mapp_tx));                                % Convert binary columns to integers, fliplr to reverse column order

        % Integer to Symbol Conversion
        sigma = sqrt(sum(([0:2^Nbps-1] - (2^Nbps-1)/2).^2) / 2^Nbps);   % Calculate the standard deviation
        symb_tx = 1/sigma * (int_tx - (2^Nbps-1)/2);                    % Normalize symbols to have unit variance

    case 'qam'

        % REAL PART
        NbpsI = Nbps / 2;                                               % Number of bits per symbol for the real part
        bit_tx2I = bit_tx2(:, 1:NbpsI);                                 % Split the matrix vertically to get the real part

        % Gray to Binary Conversion for the Real Part
        mapp_txI(:, 1) = bit_tx2I(:, 1);
        for ii = 2:NbpsI
            mapp_txI(:, ii) = xor(mapp_txI(:, ii-1), bit_tx2I(:, ii));  % XOR operation to convert Gray code to binary
        end

        % Binary to Integer Conversion for the Real Part
        int_txI = bi2de(fliplr(mapp_txI));                              % Convert binary columns to integers, fliplr to reverse column order

        % Integer to Symbol Conversion for the Real Part
        sigmaI = sqrt(sum(([0:2^NbpsI-1] - (2^NbpsI-1)/2).^2) / 2^NbpsI); % Calculate the standard deviation
        symb_txI = 1/sigmaI/sqrt(2) * (int_txI - (2^NbpsI-1)/2);          % Normalize symbols to have unit variance

        % IMAGINARY PART
        NbpsQ = Nbps / 2;                                               % Number of bits per symbol for the imaginary part
        bit_tx2Q = bit_tx2(:, NbpsQ+1:end);                             % Split the matrix vertically to get the imaginary part

        % Gray to Binary Conversion for the Imaginary Part
        mapp_txQ(:, 1) = bit_tx2Q(:, 1);
        for ii = 2:NbpsQ
            mapp_txQ(:, ii) = xor(mapp_txQ(:, ii-1), bit_tx2Q(:, ii));  % XOR operation to convert Gray code to binary
        end

        % Binary to Integer Conversion for the Imaginary Part
        int_txQ = bi2de(fliplr(mapp_txQ));                              % Convert binary columns to integers, fliplr to reverse column order

        % Integer to Symbol Conversion for the Imaginary Part
        sigmaQ = sqrt(sum(([0:2^NbpsQ-1] - (2^NbpsQ-1)/2).^2) / 2^NbpsQ); % Calculate the standard deviation
        symb_txQ = 1/sigmaQ/sqrt(2) * (int_txQ - (2^NbpsQ-1)/2);          % Normalize symbols to have unit variance

        % COMPLEX SYMBOL
        symb_tx = symb_txI + 1i * symb_txQ;                             % Combine real and imaginary parts to form complex symbols

end
