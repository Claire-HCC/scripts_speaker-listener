function [sig] = circulagregress2(x,xx,lags)
% negative shift means xx precedes, positive shift means xx follows
% quick and dirtxx function to calculate lagged correlations

Nlag = length(lags);
Nsamp = length(x);
Nxx=size(xx,2);
b = zeros(Nxx,Nlag);
X=[];
Nsamp=length(x);

for ilag = 1:Nlag
    shift = lags(ilag);

    if shift > 0
        t =[[ (1 + abs(shift) ) :  Nsamp]  [1:(abs(shift) )]  ];
        X(:,(end+1):(end+Nxx))=xx(t,:);
    elseif shift < 0
        t = [(Nsamp-abs(shift)+1):Nsamp 1:(Nsamp-abs(shift))];
        X(:,(end+1):(end+Nxx))=xx(t,:);
    elseif shift == 0
        X(:,(end+1):(end+Nxx))=xx;
    end
end

md1=stepwiselm(X,x,'Constant','Upper','Linear','criterion','rsquared');
Xi=cell2mat(cellfun(@(x) str2num(strrep(x,'x','')),md1.CoefficientNames,'UniformOutput',0));
sigi=reshape(1:Nxx*Nlag,Nxx,Nlag);
sig=zeros(size(sigi));
sig(ismember(sigi,Xi))=1;


