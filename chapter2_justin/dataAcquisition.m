function [toa, cfo] = dataAcquisition(sig, a, K, Tsymb)

    N = length(a); % pilot length
    D = zeros(length(K),length(sig)); % K is the cross-correlation averaging window

    sumK = zeros(1,length(sig));
    y = [sig zeros(1,N)]; % zero pad the signal otherwise exceed index error

    for k = 1:K
        % Compute Dk for each n
        for n = 1:length(sig)-N+1
            sumDk = 0;
            for l = k:N-1
                sumDk = sumDk + (conj(y(n+l))*a(l+1)) * conj((conj(y(n+l-k))*a(l-k+1)));
            end
            D(k,n) = (1/(N-k))*sumDk; % differential cross correlation
        end
        % Update sum to maximaze in order to find n_estimate 
        sumK = sumK + abs(D(k,:));
    end

    [~, toa] = max(sumK);
    
%     %Cpmpute CFO
%      sumDeltaF = 0;
%     for k = 1:K
%         sumDk = 0;
%             for L = k+1:N
%                 sumDk = sumDk + (conj(y(toa+L))*a(L)) * conj((conj(y(toa+L-k))*a(L-k)));
%             end
%             Dk_toa = (1/(N-k))*sumDk; %Differential cross correlation
%             sumDeltaF = sumDeltaF + angle(Dk_toa)/2*pi*k*Tsymb;
%     end
    sumDeltaF = 0;
    for k = 1:K
        sumDeltaF = sumDeltaF + ( angle(D(k,toa))/(2*pi*k*Tsymb) );
    end
    cfo = -(1/K)*sumDeltaF;

end