function [r_corrected, est_time_error] = gardner(r, K, OSF)
    % Arguments:
    %   r   : received (oversampled) signal
    %   K   : loop gain (kappa), controls timing correction step size
    %   OSF : oversampling factor (samples per symbol)

    % NOTE: Downsampling is performed during the correction process
    
    y_eps = zeros(1, length(r)/OSF);                % Preallocate corrected signal (1 sample per symbol)
    est_time_error = zeros(1, length(r)/OSF);       % Initialize epsilon
    y_eps(1) = r(1);                                % Initialize first corrected sample
    
    for n = 1:(length(r)/OSF)-1
        x_vector = ((n-1)*OSF+1 : OSF*n+2);         % Sample indices for current symbol + padding
        symbols = r(x_vector);                      % Extract current symbol window from received signal
        
        % Compute interpolation positions (fractional timing offsets)
        n_err = n*OSF + 1 - est_time_error(n);              % Target sample timing (symbol edge)
        n_mid_err = n*OSF + 1 - OSF/2 - est_time_error(n);  % Mid-symbol sample (for error computation)
        
        y_eps(n+1) = interp1(x_vector, symbols, n_err, 'linear');   % Interpolate to get y[n]
        y_mid = interp1(x_vector, symbols, n_mid_err, 'linear');    % Interpolate mid-symbol y[n-1/2]
        
        % Gardner timing error detector formula
        est_time_error(n+1) = est_time_error(n) + 2*K*real(y_mid * (conj(y_eps(n+1)) - conj(y_eps(n))));
    end

    est_time_error = est_time_error / OSF;  % Normalize time error (optional, for readability)
    r_corrected = y_eps;                    % Return corrected, symbol-rate signal
end
