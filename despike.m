function [dout] = despike(din,thr,n);
% usage: din = data in, limit is size in multiples of sigma willing to accept
% [audhrf] = despike(audhrf,4);
%clean nan's
rm = find(isnan(din));
gd = find(~isnan(din));

spike = [find(din  > thr(2)); find(din  < thr(1))];
dout = din;
dout(spike)=NaN;


% replace the bad values in data_out with the average of the neighboring
% datapoints
for spi = 1:length(spike);
    sp=spike(spi);
    dout(sp)=  nanmean(dout(max(1,min(1,sp-n/2)):max(1,min(length(dout),sp+n/2))));
    
end % for
