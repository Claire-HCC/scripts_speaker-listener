function [F,c_v] = granger_cause(y,x,alpha,lag)
% [F,c_v] = granger_cause(y,x,alpha,lags)
% Granger Causality test
% Does y Granger Cause x?
%
% User-Specified Inputs:
%   y -- A column vector of data
%   x -- A column vector of data
%   alpha -- the significance level specified by the user
%   lags -- the maximum number of lags to be considered
% User-requested Output:
%   F -- The value of the F-statistic
%   c_v -- The critical value from the F-distribution
%
% The lag length selection is chosen using the Bayesian information
% Criterion 
% Note that if F > c_v we reject the null hypothesis that x does not
% Granger Cause y
% Chandler Lutz, UCR 2009
% Questions/Comments: chandler.lutz@email.ucr.edu
% $Revision: 1.0.0 $  $Date: 09/30/2009 $
% $Revision: 1.0.1 $  $Date: 10/20/2009 $
% $Revision: 1.0.2 $  $Date: 03/18/2009 $
% References:
% [1] Granger, C.W.J., 1969. "Investigating causal relations by econometric
%     models and cross-spectral methods". Econometrica 37 (3), 424?438.
% Acknowledgements:
%   I would like to thank Mads Dyrholm for his helpful comments and
%   suggestions

%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of y
tn = size(y,2);
voxn=size(y,1);

BIC = zeros(lags,1);
%Specify a matrix for the restricted RSS
RSS_R = zeros(lags,1);

for li=1:length(lags);
    lags_temp=lags(1:li);
    
    ystar = y(:,max(abs(lags_temp))+1:T);
    xstar = [zeros(voxn,length(lags_temp))];
    %Populate the xstar matrix with the corresponding vectors of lags

    for j=1:length(lags_temp);
        lag_temp=lags_temp(j);
        xstar(:,j+1) = x(lag_temp+1-j:T-j);
    end
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [b,bint,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
    BIC(lag,:) = T*log(r'*r/T) + (lag+1)*log(T);
    
    %Put the restricted residual sum of squares in the RSS_R vector
    RSS_R(lag,:) = r'*r;
    
    lag = lag+1;
    
end
[dummy,x_lag] = min(BIC);
%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of x
BIC = zeros(lags,1);
%Specify a matrix for the unrestricted RSS
RSS_U = zeros(lags,1);
lag = 1;
while lag <= lags
    
    ystar = y(lag+x_lag+1:T,:);
    xstar = [ones(T-(lag+x_lag),1) zeros(T-(lag+x_lag),x_lag+lag)];
    %Populate the xstar matrix with the corresponding vectors of lags of y
    j = 1;
    while j <= x_lag
        xstar(:,j+1) = y(lag+x_lag+1-j:T-j,:);
        j = j+1;
    end
    %Populate the xstar matrix with the corresponding vectors of lags of x
    j = 1;
    while j <= lag
        xstar(:,x_lag+j+1) = x(lag+x_lag+1-j:T-j,:);
        j = j+1;
    end
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [b,bint,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
    BIC(lag,:) = T*log(r'*r/T) + (lag+1)*log(T);
    
    RSS_U(lag,:) = r'*r;
    
    lag = lag+1;
    
end
[dummy,y_lag] =min(BIC);
%The numerator of the F-statistic
F_num = ((RSS_R(x_lag,:) - RSS_U(y_lag,:))/y_lag);
%The denominator of the F-statistic
F_den = RSS_U(y_lag,:)/(T-(x_lag+y_lag+1));
%The F-Statistic
F = F_num/F_den;
c_v = finv(1-alpha,y_lag,(T-(x_lag+y_lag+1)));
    
    
    
    
