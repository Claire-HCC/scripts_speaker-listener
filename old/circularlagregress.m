function [b] = circulagregress(x,xx,lags)
% negative shift means xx precedes, positive shift means xx follows
% quick and dirtxx function to calculate lagged correlations

Nlag = length(lags);
Nsamp = length(x);

b = zeros(size(xx,2),Nlag);

Nsamp=length(x);

for ilag = 1:Nlag
    shift = lags(ilag);
    
    t1 = [1 : Nsamp  ];
    t2 =[[ (1 + abs(shift) ) :  Nsamp]  [1:(abs(shift) )]  ];
    
    if shift > 0
        b(:,ilag) = regress(x(t1), xx(t2,:));
    elseif shift < 0
        b(:,ilag) = regress(x(t2), xx(t1,:));
    elseif shift == 0
        b(:,ilag) = regress(x,xx);
    else
        b(:,ilag) = NaN;
    end
    
end

    
    
    
