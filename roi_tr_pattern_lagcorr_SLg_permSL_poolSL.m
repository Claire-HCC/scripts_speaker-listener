
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

lags_tested={-10:10,  -30:30};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        r_perm=nan([length(rnames)  length(lags) ]);
        for perm=1:48;
            
            
            if exist(sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/SLg/perm/SL_lag%d-%dpermSL%02d.mat' , expdir,exp, timeUnit,froidir,min(lags),max(lags),perm));
                load( sprintf('%s/%s/fmri/pattern_lagcorr/%s/roi/%s/SLg/perm/SL_lag%d-%dpermSL%02d' , expdir,exp, timeUnit,froidir,min(lags),max(lags),perm),'r','lags','rnames','keptT');
                r_perm(:,:,perm)=r;
            end
            
            r=r_perm;
            save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permSL' ],'rnames','r','lags','keptT');
        end
    end
end

