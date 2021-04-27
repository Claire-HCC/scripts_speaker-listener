
% use cluster from pieman old to replace mor.
% i probably will have to define my own roi...
loc='mypc';
set_parameters;
timeUnit='tr' ;
exp='pieman_rest';
ei=find(ismember(exp_parameters.experiments,exp));
dice=[];

for k=5;
    perc=0.3;
    kmeansd='cityblock';
    
    load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_g2'],'rzm_g1','rzm_g2','keptvox','keptT');
    
    [idx_g1]= kmeans(rzm_g1,k,'Distance',kmeansd);
    [idx_g2]= kmeans(rzm_g2,k,'Distance',kmeansd);
    
    
    for ki_g1=1:k;
        for ki_g2=1:k;
            dice(ki_g1,ki_g2,k)=2*sum(idx_g1==ki_g1 & idx_g2==ki_g2)/(sum(idx_g1==ki_g1) + sum(idx_g2==ki_g2));
        end
    end
end

