%% Power Density Spectrum
%
% Author: Tilemachos S. Doganis
%
% To be run after fsk_4.m and psk_4.m
%
% Visualizes the Power Density Spectrum of the 4-FSK and 4-PSK encoded
% signals. Due to the fact that more bandwidth is used for FSK, power is
% more widely distributed near the half of the carrier frequency (51.2 kHz)
% while in PSK there is a sharper power drop. Considering the cost in
% bandwitdh, PSK seems more cost effective.
%
T_sym = 4*10^(-6);
N = 2048;
Fs = N/(T_sym*(L_b/2));
freq = linspace(-Fs,Fs,numel(Ps_m));

PSD = 1/length(Ps_m) * abs(fft(Ps_m)).^2;
FSD = 1/length(Fs_m) * abs(fft(Fs_m)).^2;

figure, 
semilogy(freq,PSD), grid on, hold all, title('Power Density Spectrum')
semilogy(freq,FSD)
xlabel('Frequency (Hz)') % x-axis label
ylabel('Power Density (W/Hz)') % y-axis label
legend('PSK','FSK')
xlim([-Fs Fs])