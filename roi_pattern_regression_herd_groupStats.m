clear all
set_parameters;

timeUnit='tr';
froidir='mor';
lags=-4:4;
for ei=1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_'  exp '_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herd','rnames');
    herds_z(:,ei)=0.5*log((1+herd)./(1-herd));
end

[~,p,~,stats]=ttest(herds_z');
