function [lagcc, lagccblockvals] = circularlagcorr_byCol(x,y,lags)
% positive shift means y follows, negative shift means y precedes
% quick and dirty function to calculate lagged correlations
% lagcc = values of correlation coefficients
if nargin < 3; lags = -10:10; end


Nlag = length(lags);
Nsamp = size(x,1);

lagcc = zeros(Nlag,size(x,2));

for ilag = 1:Nlag
    shift = lags(ilag);
    
    t1 = [1 : Nsamp  ];
    t2 =[[ (1 + abs(shift) ) :  Nsamp]  [1:(abs(shift) )]  ];
    
    if shift > 0
        lagcc(ilag,:) = corr_col(x(t1,:), y(t2,:));
    elseif shift < 0
        lagcc(ilag,:) = corr_col(x(t2,:), y(t1,:));
    elseif shift == 0
        lagcc(ilag,:) = corr_col(x,y);
    else
        lagcc = NaN;
    end
end


