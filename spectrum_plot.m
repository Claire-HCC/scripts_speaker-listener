function spectrum_plot(x,fs);
% 
% m=length(x);
% n = pow2(nextpow2(m));  % transform length
% y = fft(x,n);        % DFT of signal
% f = (0:n-1)*(fs/n);
% power = abs(y).^2/n;      
% 
% figure
% plot(f(1:floor(n/2)),power(1:floor(n/2)))
% xlabel('Frequency')
% ylabel('Power')


n = length(x);
xdft = fft(x);
xdft = xdft(1:n/2+1);
psdx = (1/(fs*n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
end