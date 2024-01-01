function smoothed_data = KalmanSmoother(input_data)
    % 
    A = 1;  % 
    H = 1;  % 
    Q = 1e-4;  % 
    R = 1e-2;  % 

    % 
    x_hat = input_data(1);  % 
    P = 1;  % 

    % 
    for k = 2:length(input_data)
        % 
        x_hat_minus = A * x_hat;
        P_minus = A * P * A' + Q;

        % 
        K = P_minus * H' / (H * P_minus * H' + R);
        x_hat = x_hat_minus + K * (input_data(k) - H * x_hat_minus);
        P = (1 - K * H) * P_minus;

        % 
        smoothed_data(k) = x_hat;
    end
end
