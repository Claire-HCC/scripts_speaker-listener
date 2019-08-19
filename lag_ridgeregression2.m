function [ B_regAll, B_regAuto, rss_regAll, rss_regAuto] = lag_ridgeregression2(speaker, listeners, window, RR_exp);
% ridfe regression with the past datapoint in y to solve the
% auto-correlation problem
% win must be < 0
% This is similar granger causality
% Computes coupling between speaker voxel timecourses and timecourses from
% another listeners using a linear filter approach. For example, to compute the
% coupling between a speaker and a listeners of listeners, the speaker would be
% the voxel timecourses from the speaker and the listeners would be the voxel
% timecourses from all listeners. The 2-element window parameter specifies
% the starting and ending offsets for the speaker. RR_exp is an optional
% parameter that adds ridge regularization to the beta estimates, which
% empirically does not appear necessary for windows that are relatively
% small compared to the length of the timecourse.

% It can return the average fitted coupling filters in B. The absolute scale of
% these values is uninformative, but the sign and relative scale is
% meaningful.

% speaker = voxels by time
% listeners = voxels by time by subject
% window = [min max] speaker offset
% RR_exp (optional) sets ridge regression strength to 10^RR_exp

% coupling = voxels by subject
% B_regAll = ridge coefficients. voxels by window size. This is obtained by
% including all time lagged time series in one model.
% B_regOne = ridge coefficients. voxels by window size. This is obtained by
% construct a one regresor model for each lag.

if (nargin < 4)
    RR_exp = -Inf;
end

nv = size(speaker,1);
ns = size(listeners,3);
T = size(speaker,2);

listeners_time = window(1):window(2);
B = zeros(nv,length(listeners_time));
% nsamp = T-length(listeners_time)+1;

minT = 1+max(-1*window(1),0);
maxT = T-max(window(2),0);
nsamp = (maxT-minT+1);

for v = 1:nv
    y = zscore(speaker(v,minT:maxT))';
    
    Y = zeros(nsamp,length(listeners_time));
    for t = minT:maxT
        Y(t-minT+1,:) = speaker(v,t+listeners_time);
    end
    Y = zscore(Y); % X: time series x lags
     Y(:,listeners_time==0)=[];
     
    for loo_subj = 1:ns;
        X = zeros(nsamp,length(listeners_time));
        for t = minT:maxT
            X(t-minT+1,:) = listeners(v,t+listeners_time,loo_subj);
        end
        X = zscore(X); % X: time series x lags
        X(:,listeners_time==0)=[];
        
        cv_B = ridge(y,X,10^RR_exp);
        rss_regAll(v,1,loo_subj)=(y-X*cv_B)'*(y-X*cv_B);
        B_regAll(v,:,loo_subj) = zscore(cv_B); % zscore
        
        cv_B = ridge(y,Y,10^RR_exp);
        rss_regAuto(v,1,loo_subj)=(y-Y*cv_B)'*(y-Y*cv_B);
        B_regAuto(v,:,loo_subj) = zscore(cv_B); % zscore
        
    end
    
end
