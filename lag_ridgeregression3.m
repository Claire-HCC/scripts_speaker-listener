function [ B, F, p] = lag_ridgeregression3(dv, iv, window, RR_exp);
% ridge regression with the future datapoint in iv to explain dv
% dv = voxels by time
% iv = voxels by time by subject
% window = [min max]
% RR_exp (optional) sets ridge regression strength to 10^RR_exp


if (nargin < 4)
    RR_exp = -Inf;
end

nv = size(dv,1);
T = size(dv,2);

iv_time = window(1):window(2);
B = zeros(nv,length(iv_time));
% nsamp = T-length(iv_time)+1;

minT = 1+max(-1*window(1),0);
maxT = T-max(window(2),0);
nsamp = (maxT-minT+1);

for v = 1:nv;
    
    y = zscore(dv(v,minT:maxT))';
    
    X = zeros(nsamp,length(iv_time));
    for t = minT:maxT
        X(t-minT+1,:) = iv(v,t+iv_time);
    end
    X = zscore(X); % X: time series x lags
    
    cv_B = ridge(y,X,10^RR_exp);
    dfm=length(cv_B);
    dfe=(length(y)-1-length(cv_B));
    msm=(X*cv_B-mean(y))'*(X*cv_B-mean(y))/dfm;
    mse=(y-X*cv_B)'*(y-X*cv_B)/dfe;
    F(v,1)=msm/mse;
    p(v,1)=1-fcdf(F(v,1),dfm,dfe);
    
    B(v,:) = zscore((cv_B)); % zscore makes a huage difference
   
    
end


