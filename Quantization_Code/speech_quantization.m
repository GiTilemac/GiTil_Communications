%% Voice Signal Quantization
%
% Author: Tilemachos S. Doganis
%
% Load a .wav audio file and quantize it using Uniform PCM and Lloyd-Max
% with 2, 4 and 6 bits/sample
%
close all
clear
clc

[y,fs] = audioread('speech.wav');
N = 2;
fprintf('Quantizing voice signal using Lloyd-Max and %d bits...\n',N);
[yq2, centers2, D2, SQNR2] = Lloyd_Max(y, N, -1, 1);

N = 4;
fprintf('Quantizing voice signal using Lloyd-Max and %d bits...\n',N);
[yq4, centers4, D4, SQNR4] = Lloyd_Max(y, N, -1, 1);

N = 6;
fprintf('Quantizing voice signal using Lloyd-Max and %d bits...\n',N);
[yq6, centers6, D6, SQNR6] = Lloyd_Max(y, N, -1, 1);

figure(99), hold on, grid on
plot(SQNR2,'b')
plot(SQNR4,'g')
plot(SQNR6,'r')
legend('N = 2','N = 4', 'N = 6','location','Southeast')
title('SQNR Curves')
ylabel('SQNR (dB)')
xlabel('Cycles')