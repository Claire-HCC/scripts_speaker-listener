function [lagcc] = circularlagcorr_spatialTemporal_canonical(x,y,lags)
% lagcorr between two matrix by column
% positive lag means y follows
% negative lag means x preceds

% quick and dirty function to calculate lagged correlations
% lagcc = values of correlation coefficients
%
% Nblocks = number of blocks into which to divide data, if you want to
% calculate correlations in blocks and take the median. default = 1 (no
% subdivision of the data )


if nargin < 3; lags = -10:10; end

Nlag = length(lags);
Nsamp = size(y,2);

lagcc = zeros(1,Nlag);

for ilag = 1:Nlag
    shift = lags(ilag);
    t1 = [1 : Nsamp  ];
    t2 =[[ (1 + abs(shift) ) :  Nsamp]  [1:(abs(shift) )]  ];
    
    if shift > 0
        y_temp=y(:,t1);
        x_temp=x(:,t2);
        [~,~, r] = canoncorr(x_temp',y_temp');
        lagcc(1,ilag) =  r(1);
    elseif shift < 0
        y_temp=y(:,t2);
        x_temp=x(:,t1);
        [~,~, r] = canoncorr(x_temp',y_temp');
        lagcc(1,ilag) =  r(1);
    elseif shift == 0
        [~,~, r] = canoncorr(x',y');
        lagcc(1,ilag) =  r(1);
    else
        lagcc = NaN;
    end
    
    
end


