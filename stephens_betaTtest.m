
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-4:4;
contrast=[ ]

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','keptvox','p');
        pfdr=nan(size(p));
        [~,~,pfdr]=fdr(p);
        
        
y=rand(100,1);
X=rand(100,10);
X1=[ones(100,1) X];
c=[0 1 zeros(1,9)]';
[b,bint,r,rint,stats] = regress(y,[X1]);
md1=fitlm(X,y);
md1
sum(c.*b)/sqrt(stats(4)*(c'*inv(X1'*X1)*c))

% stats(4) is error variance, which can be computed from residuals:  sum(r.^2)/(length(y)-size(X1,2))