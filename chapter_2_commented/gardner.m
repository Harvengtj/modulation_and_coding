function [r_corrected,est_time_error] = gardner(r,K,OSF)

    % Note: Downsampling is performed implicitly within this loop
    
    y_eps = zeros(1,length(r)/OSF);                  % Output signal (timing corrected)
    est_time_error = zeros(1,length(r)/OSF);         % Estimated timing error per symbol
    y_eps(1)=r(1);                                    % Initialize first output sample
    
    for n = 1:length(r)/OSF-1
        x_vector=((n-1)*OSF+1:OSF*n+2);               % Sample indices for interpolation window
        symbols = r(x_vector);                        % Extract window of samples
        n_err = n*OSF+1-est_time_error(n);            % Time index for current symbol
        n_mid_err = n*OSF+1-OSF/2 -est_time_error(n); % Time index for midpoint between symbols
        
        y_eps(n+1) = interp1(x_vector,symbols,n_err,'linear');     % Interpolated symbol y[n]
        y_mid = interp1(x_vector,symbols,n_mid_err,'linear');      % Interpolated midpoint y[n - 1/2]
        
        % Update estimated timing error using Gardner's algorithm
        est_time_error(n+1)=est_time_error(n)+2*K*real(y_mid*(conj(y_eps(n+1))-conj(y_eps(n))));
    end

    est_time_error=est_time_error/OSF;                % Normalize error by OSF
    r_corrected = y_eps;                              % Return corrected signal
end
