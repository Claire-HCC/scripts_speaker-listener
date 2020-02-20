function [lagcc] = lagcorr_spatialTemporal(x,y,lags)
% lagcorr between two matrix by column
% positive lag means x precedes, y follows (lag=1, corr(x(1:(end-1)),y(2:end))).
% negative lag means series x follows, y preceds (lag=-1, corr(x(2:(end)),y(1:(end-1)))).

% quick and dirty function to calculate lagged correlations
% lagcc = values of correlation coefficients
%
% Nblocks = number of blocks into which to divide data, if you want to
% calculate correlations in blocks and take the median. default = 1 (no
% subdivision of the data )


if nargin < 3; lags = -10:10; end

Nlag = length(lags);
Nsamp = size(x,2);

lagcc = zeros(1,Nlag);

for ilag = 1:Nlag
    shift = lags(ilag);
    
    t1 = 1 : (Nsamp - abs(shift));
    t2 = (1 + abs(shift) ) :  Nsamp;
    
    if shift > 0
        x_temp=x(:,t1);
        y_temp=y(:,t2);
        lagcc(1,ilag) = corr(x_temp(:),y_temp(:));
    elseif shift < 0
       x_temp=x(:,t2);
        y_temp=y(:,t1);
        lagcc(1,ilag) = corr(x_temp(:),y_temp(:));
    elseif shift == 0
        lagcc(1,ilag) = corr(x(:),y(:));
    else
        lagcc = NaN;
    end
    
    
end


