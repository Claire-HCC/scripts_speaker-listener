function [f power]=spectrum_plot(x,fs);


% plot waveform
t=(0:length(x)-1)/fs;        % times of sampling instants
subplot(4,1,1);
plot(t,x);
legend('Waveform');
xlabel('Time (s)');
ylabel('Amplitude');


m=length(x);
n = pow2(nextpow2(m));  % transform length
y = fft(x,n);        % DFT of signal
f = (0:n-1)*(fs/n);
power = abs(y).^2/n;      

subplot(4,1,2);
f=f(1:floor(n/2));
power=power(1:floor(n/2));
plot(f,power)
xlabel('Frequency')
ylabel('Power')


n = length(x);
xdft = fft(x);
xdft = xdft(1:n/2+1);
psdx = (1/(fs*n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

subplot(4,1,3);
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');


% do fourier transform of windowed signal
Y=xdft;
% cepstrum is DFT of log spectrum
C=fft(log(abs(Y)+eps));
% plot between 1ms (=1000Hz) and 20ms (=50Hz)
%ms1=fs/1;                 % maximum speech Fx at 50Hz
% ms20=fs/0.01;                  % minimum speech Fx at 1Hz
%q=(ms1:ms20)/fs;
q=(1:length(C))/fs;
%plot(q,abs(C(ms1:ms20)));
subplot(4,1,4);
plot(q,abs(C));
 xlim([0.5 100]);
ylim([0 50]);
legend('Cepstrum');
xlabel('Quefrency (s)');
ylabel('Amplitude');

end

