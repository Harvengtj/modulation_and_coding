function [toa, cfo] = dataAcquisition(sig, a, K, Tsymb)

    N = length(a);                                  % Length of known preamble sequence
    D = zeros(length(K),length(sig));               % Cross-correlation matrix

    sumK = zeros(1,length(sig));                    % Accumulator for correlation magnitudes
    y = [sig zeros(1,N)];                           % Zero-pad signal to prevent indexing errors

    for k = 1:K
        % Compute differential correlation D_k[n] for each time shift n
        for n = 1:length(sig)-N+1
            sumDk = 0;
            for L = k:N-1
                % Differential correlation using delayed versions of signal and known preamble
                sumDk = sumDk + (conj(y(n+L))*a(L+1)) * conj((conj(y(n+L-k))*a(L-k+1)));
            end
            D(k,n) = (1/(N-k)) * sumDk;             % Normalize correlation
        end
        % Accumulate magnitude of D_k[n] across different k
        sumK = sumK + abs(D(k,:));
    end

    [~, toa] = max(sumK);                           % Time of arrival estimated by peak correlation
    
    sumDeltaF = 0;
    for k = 1:K
        % Estimate frequency offset using phase of correlation
        sumDeltaF = sumDeltaF + ( angle(D(k,toa)) / (2*pi*k*Tsymb) );
    end
    cfo = -(1/K) * sumDeltaF;                       % Average to estimate carrier frequency offset

end
