function [symb_tx] = mapping(bit_tx,Nbps,modulation)

% INPUTS:
% - bit_tx : vector of input bits 
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - symb_tx : vector of ouput symbols (variance 1)


%% STEP 1
Nsymb = size(bit_tx,1)/Nbps; % Number of symbols

bit_tx2 = reshape(bit_tx,Nbps,Nsymb)';

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
switch modulation,
    
    case 'pam'
        
        % Gray to binary
        mapp_tx(:,1) = bit_tx2(:,1);
        for ii = 2:Nbps,
           mapp_tx(:,ii) = xor( mapp_tx(:,ii-1) , bit_tx2(:,ii) ); 
        end

        % Binary to integer
        int_tx = bi2de(fliplr(mapp_tx)); % fliplr --> invert columns order to put the most significant bits to the left
                                         % bi2de --> converts binary bits to integer : [3;2;1;0]

        % Integer to symbol
        sigma = sqrt(sum(([0:2^Nbps-1]-(2^Nbps-1)/2).^2)/2^Nbps); % sqrt[E[(X-E[X])]^2] = 0.79
        symb_tx = 1/sigma * (int_tx - (2^Nbps-1)/2); % [1.9; 0.63; -0.63; -1.9] (normalized symbols)
        
    case 'qam'
        
        % REAL PART
        NbpsI = Nbps/2;
        bit_tx2I = bit_tx2(:,1:NbpsI); % Cut vertically the matrix in the middle to separate real and imaginary parts
        
        % Gray to binary
        mapp_txI(:,1) = bit_tx2I(:,1);
        for ii = 2:NbpsI,
           mapp_txI(:,ii) = xor( mapp_txI(:,ii-1) , bit_tx2I(:,ii) ); 
        end

        % Binary to integer
        int_txI = bi2de(fliplr(mapp_txI));

        % Integer to symbol
        sigmaI = sqrt(sum(([0:2^NbpsI-1]-(2^NbpsI-1)/2).^2)/2^NbpsI); 
        symb_txI = 1/sigmaI/sqrt(2) * (int_txI - (2^NbpsI-1)/2);
        
        
        % IMAGINARY PART
        NbpsQ = Nbps/2;
        bit_tx2Q = bit_tx2(:,NbpsQ+1:end);
        
        % Gray to binary
        mapp_txQ(:,1) = bit_tx2Q(:,1);
        for ii = 2:NbpsQ,
           mapp_txQ(:,ii) = xor( mapp_txQ(:,ii-1) , bit_tx2Q(:,ii) ); 
        end

        % Binary to integer
        int_txQ = bi2de(fliplr(mapp_txQ));

        % Integer to symbol
        sigmaQ = sqrt(sum(([0:2^NbpsQ-1]-(2^NbpsQ-1)/2).^2)/2^NbpsQ); 
        symb_txQ = 1/sigmaQ/sqrt(2) * (int_txQ - (2^NbpsQ-1)/2);
       
        
        % COMPLEX SYMBOL
        symb_tx = symb_txI + j*symb_txQ;
       
end