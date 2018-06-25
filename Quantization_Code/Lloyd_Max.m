%% Lloyd-Max Quantization Function
%
% Author: Tilemachos S. Doganis
%
% Input: Signal 'x', bits/sample 'N', lower signal amplitude bound 'min_value',
% upper signal bound 'max_value'
% Output: quantization centers vector 'centers',
% sample - level matching indices 'xq',
% mean distortion per step vector 'D',
% Signal to Quantization Noise ratio 'SQNR'
%
%   Initially, uniform quantization is performed for comparison purposes.
% Quantization region boundaries 'T' are calculated as the center means,
% SQNR and distortion are stored in the first position of 'SQNR' and 'D',
% respectively.
%
%   The Lloyd-Max loop begins now. On each step, the difference between the  
% last two consecutive distortions is checked, and if it is smaller than
% threshold 'eps', the algorithm has converged and stops.
% Afterwards, the arithmetic mean of each region's samples is calculated
% and is designated as the new center, and the new region bounds are placed
% in the middle between the aforementioned centers. 
% Lastly, the assigned quantization levels are stored in 'xq' and the
% distortion difference between the current and previous levels is stored
% in 'D'. 
% 
function [xq, centers, D, SQNR] = Lloyd_Max(x, N, min_value, max_value)
%% Initialization
n = numel(x);
m = 2^N;

% Quantization Region Bounds T
T = zeros(m-1,1);
count = zeros(m,1);

% Sort sample values
[y, I] = sortrows(x);

% Quantization Region Size Delta
rsize = abs(max_value - min_value)/m;
a = min_value;

%% Bound values between min_value and max_value
% Bound minimum values
i = 1;
while y(i) < min_value
    y(i) = min_value;
    x(I(i)) = min_value;
    i = i+1;
end

% Bound maximum values
i = n;
while y(i) > max_value
    y(i) = max_value;
    x(I(i)) = max_value;
    i = i-1;
end

%% Uniform PCM
%Calculate uniform region centers
centers = zeros(m,1);
for i = 1:m
   centers(i) = (2*a+(2*i-1)*rsize)/2;
end

%Calculate bounds T
for k = 1:m-1
   T(k) = ( centers(k) + centers(k+1) )/2;
end

%% Uniform quantization (for comparison)
%Store uniformaly quantized region for each x
xq = zeros(n,1);
i = 1;
flag = 0;
for k = 1:m-1
    while ((y(i) < T(k)) && (flag == 0))
        xq(I(i))= k;
        if i < n
            i = i+1;
        else
            flag = 1;
        end
    end
end
xq(I(i:n)) = m;

%Plot Voice Signal and its Uniform Quantization
figure
subplot(3,1,1), plot(1:numel(x),x,'r');
title('Voice Signal');
subplot(3,1,2), plot(1:numel(x),centers(xq),'b');
title(['Uniform PCM (N = ' num2str(N) ' bits/sample)']);

%Calculate Distortion
D(1) = sum((y-centers(xq(I))).^2)/n;

%Calculate Uniform SQNR
SQNR(1) = 10*log10(sum(y.^2)/sum((y-centers(xq(I))).^2));


%% Lloyd-Max Algorithm
dfr = 1;
Kmax = 1;
while dfr >= eps
    %Calculate new centers
    i = 1;
    flag = 0;
    
    %Find mean value of region k
    for k = 1:m-1
        count(k) = 0;
        s = 0;
        while (y(i) < T(k) && flag == 0)
            count(k) = count(k) + 1;
            s = s+y(i);
            if (i>n-1)
                flag = 1;
            else
                i = i+1;
            end
        end
        %New center is mean value of region
        if ( count(k) > 0 )
            centers(k) = s / count(k);
        end   
    end
    %Find last center
    if (flag == 0)
        count(m) = n-i+1;
        s = sum(y(i:n));
        centers(m) = s / count(m);
    end

    %Calculate new bounds
    for k=1:m-1
       T(k) = ( centers(k) + centers(k+1) )/2;
    end

    %Store each Quantized Region j to corresponding index xq
    i = 1;
    flag = 0;
    for k = 1:m-1
        while ((y(i) < T(k)) && (flag == 0))
            xq(I(i))= k;
            if (i<n)
                i = i+1;
            else
                flag = 1;
            end
        end
    end
    xq(I(i:n))= m;

    %Calculate Distortion
    D(Kmax+1) = sum((y-centers(xq(I))).^2)/n;
    
    %Calculate SQNR
    SQNR(Kmax+1) = 10*log10(sum(y.^2)/sum((y-centers(xq(I))).^2));
    
    dfr = abs(D(Kmax) - D(Kmax+1));
    Kmax = Kmax+1;
end

% Calculate Entropy
count = count/numel(x);

% Zero Probability Workaround
count(count == 0) = 1; %P = 1, H = 0
H = -sum(count.*log2(count));

% Theoretical Distortion Limit
D_lim = var(x)*2^(-2*N);

% Plot Non-Uniformaly Quantized Signal
subplot(3,1,3), plot(1:numel(x),centers(xq),'g')
title(['Lloyd-Max (N = ' num2str(N) ' bits/sample)']);

fprintf('Success after %d loops\n',Kmax);
fprintf('Quantizer Entropy: %f\n', H);
fprintf('MSE Distortion (Lloyd-Max): %f\n', D(end));
fprintf('MSE Distortion (Uniform PCM): %f\n',D(1));
fprintf('Distortion Limit: %f\n',D_lim);
fprintf('SQNR (Lloyd-Max): %f dB\n',SQNR(end));
fprintf('SQNR (Uniform): %f dB\n\n',SQNR(1));
end