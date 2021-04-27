
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -40:40};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        r_perm=nan([length(rnames)  length(lags) ]);
        for perm=1:48;
            
            if exist(sprintf('%s/%s/fmri/pattern/lagcorr/%s/roi/%s/SL_g/perm/lag%d-%d_permSL%02d.mat' , expdir,exp, timeUnit,froidir,min(lags),max(lags),perm));
                load( sprintf('%s/%s/fmri/pattern/lagcorr/%s/roi/%s/SL_g/perm/lag%d-%d_permSL%02d' , expdir,exp, timeUnit,froidir,min(lags),max(lags),perm),'r','lags','rnames','keptT');
                r_perm(:,:,perm)=r;
            end
            
            r=r_perm;
            save([expdir '/' exp '/fmri/pattern/lagcorr/' timeUnit '/roi/' froidir '/SL_g/perm/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL' ],'rnames','r','lags','keptT');
        end
    end
end

