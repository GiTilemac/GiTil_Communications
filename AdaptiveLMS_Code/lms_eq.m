%% Adaptive Filter Equalizer using LMS algorithm
%
% Author: Tilemachos S. Doganis
%
% INPUT
% y: Channel output after AWGN
% a: Transmitted symbol sequence
% A_m: 4-PAM dictionary
% Q: Training sequence length
% P: Total packet length
% L: Channel length
% c: Equalizer vector
% d: Learning step parameter
%
% OUTPUT
% MSE: Mean Square Root Error curve
% z_Q: Equalizer output
% c: Updated equalizer vector
%
%% Program description
% Phase 1: Training sequence
% The first 'Q' samples sent are used for training the adaptive equalizer.
% For each training symbol:
% 1) Calculate the indices of the symbols involved in ISI 'ky'
% 2) Calculate the equalization coefficients 'kc'
% 3) Implement equalization by multiplying 'kc' with 'ky'.
% 4) Calculate momentary error by removing the approximated value from the
% known one
% 5) Update the equalizer coefficients using the calculated error.
%
% Phase 2: Test sequence
% The process is the same, but because the true symbol values aren't known
% a priori, a Maximum Likelihood Detector is used, comparing the Euclidean
% norm between the equalized value 'z' and the symbol amplitude dictionary
% 'A_m'.
% The momentary error is the difference between the received and the
% equalized value.
% 

function [ MSE, zQ, c ] = lms_eq( y, a, Am, Q, P, L, c, d)
    %% Initialization
    
    % Channel length
    N = floor(L/2);

    % Output MSE
    MSE = zeros(1,P);

    % Initialize estimation
    z = zeros(P,1);
    zQ = zeros(P-Q,1);

    % Estimate power of equalizer input and calculate step Delta
    Py = (1/P)*sum(abs(y).^2);
    Delta = 1/(5*L*Py);
    Delta = d * Delta; % Optional change in proposed step

    %% LMS Training
    for k = 1:Q
        
        % Effective sample range
        ky = max(k-N,1):(k+N);
        kc = min(L,k+N):-1:1;  % 11<-6:1 

        % Filter received sample
        z(k) = c(kc)'*y(ky);

        % Calculate Error
        e = a(k) - z(k);

        % Calculate MSE
        if k > 1
            MSE(k) = (1-1/k)*MSE(k-1) + (1/k)*abs(e)^2;
        else
            MSE(k) = abs(e)^2;
        end
        
        % Update Equalizer Coefficients
        c(kc) = c(kc) + Delta*y(ky)*e';

    end

    %% LMS Core
    for k = Q+1:P
        
        % Effective sample range
        ky = (k-N):min(k+N,P);
        kc = L:-1:max(1,k-P+N+1);
        
        % Filter received sample
        z(k) = c(kc)'*y(ky);
        
        % Detector
        [~,dcs] = min(abs(Am - z(k)));
        zQ(k-Q) = Am(dcs);
        
        % Calculate Error
        e = zQ(k-Q) - z(k);
        
        % Calculate MSE
        MSE(k) = (1-1/k)*MSE(k-1) + (1/k)*abs(e)^2;

        % Update Equalizer Coefficients
        c(kc) = c(kc) + Delta*y(ky)*e';
    end
end
