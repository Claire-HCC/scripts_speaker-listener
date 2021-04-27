% Within-Cluster-Sum of Squared Errors keeps decreasing.
% wss might be good for defining the number of roi but not at the network
% level.


loc='mypc';
set_parameters;
timeUnit='tr' ;
exp='pieman_rest';
ei=find(ismember(exp_parameters.experiments,exp));
k=5;
perc=0.15;
kmeansd='cityblock';

load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked'],'rzm','keptvox','keptT');
wss=nan(1,20);
for k=11:20;
    [idx, C,sumd, D]= kmeans(rzm,k,'Distance',kmeansd);
    wss(k)=sum(sumd);
end
save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_optimalK' ],'wss','-v7.3');
