%% Random Signal Uniform Quantization
%
% Author: Tilemachos S. Doganis
%
% Produce a random complex signal with unit variance, then apply uniform
% quantization to its amplitude using 4 and then 6 bits, bounded between 0
% and 4.
%
clear
close all
clc
M = 10000;      % Signal length (in samples)
N = 4;          % Bits available for quantization
t = (randn(M,1)+1i*randn(M,1))/sqrt(2); % Produce random signal
x = abs(t).^2;  % Signal amplitude

% Quantize with N = 4 bits
fprintf('Quantizing %d samples using %d bits...\n',M,N);
[xq4, centers4] = my_quantizer(x, N, 0, 4);

% Quantize with N = 6 bits
N = 6;
fprintf('Quantizing %d samples using %d bits...\n',M,N);
[xq6, centers6] = my_quantizer(x, N, 0, 4);

figure(99), hold on, grid on
title('Sorted Samples')
plot(sort(x))
plot(sort(centers4(xq4)),'g')
plot(sort(centers6(xq6)),'r')
legend('Original Values','Quantization N = 4','Quantization N = 6')