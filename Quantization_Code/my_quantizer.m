%% Uniform Quantization Function
%
% Author: Tilemachos S. Doganis
%
function [ xq, centers ] = my_quantizer(x, N, min_value, max_value )
%
% Input: Input signal, number of bits/sample, minimum amplitude, maximum
% amplitude
% Output: indices matching samples to levels, Quantization level centers
%
% Begin by storing the sorted signal values in a vector 'y'.
% Bound the signal between a user-defined maximum and minimum values and
% count the outliers to later calculate overload probability.
% As the samples of 'y' are sorted in ascending order, assign the
% quantization center value to all samples belonging to the
% corresponding region. The permutation vector 'I' is used to assign
% maintain the original sample order.
% 
%% Initialization
n = numel(x);         % Input signal length
[y, I] = sortrows(x); % Sort samples of x in ascending order

%% Set lower and upper boundaries for sample amplitudes
% Count outliers to estimate overload probability
outliers = 0;

% Lower Bound
outliers = outliers + sum(y < min_value);
y(y < min_value) = min_value;

% Upper Bound
outliers = outliers + sum(y > max_value);
y(y > max_value) = max_value;

%% Quantization
% Calculate Region Size (Delta)
m = 2^N;                              % Number of quantization regions possible using N bits
rsize = abs(max_value - min_value)/m; % Split amplitude range into regions of same size
a = min_value;

% Calculate region centers
centers = zeros(m,1);
for i = 1:m
   centers(i) = (2*a+(i*2-1)*rsize)/2;
end

% Assign quantization values to samples
xq = zeros(n,1);
i = 1;
flag = 0;
for j = 1:m-1
    while (y(i) < a+j*rsize) && (flag == 0)
        xq(I(i)) = j;

        if i < n
            i = i+1;
        else
            flag = 1;
        end
    end
end
xq(I(i:n))= m; % Final values belong to same, final level (because y sorted)

%% Output information

% Plot Signals
figure,
subplot(2,1,1), plot(x), title('Original Signal')
subplot(2,1,2), plot(xq), title(['Quantized Signal (N = ' num2str(N) ' bits/sample)'])

% Calculate Distortion Overload Probability
if outliers > 0
    fprintf('Distortion Overload Probability: %.2f%%\n',100*outliers/n);
else
    fprintf('Distortion Overload Probability: 0%%\n'); 
end

% Experimental MSE Distortion
D = sum((y-centers(xq(I))).^2)/n;

% SQNR
SQNR = 10*log10(sum(y.^2)/sum((y-centers(xq(I))).^2));

% Theoretical Distortion
Dt = rsize^2/12;

% Theoretical Distortion Limit
D_lim = 2^(-2*N);

fprintf('SQNR: %f dB\n',SQNR);
fprintf('MSE Distortion: %f\n',D);
fprintf('Theoretical Distortion: %f\n', Dt);
fprintf('Distortion Limit: %f\n\n',D_lim);
end