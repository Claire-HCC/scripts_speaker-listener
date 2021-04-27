

% loc='mypc';
set_parameters;
timeUnit='tr' ;

froidir='isc_peak';
rois={'HG_L','precuneus','Angular_L'};


tic % 15 min
for ei=1:11;
    exp=exp_parameters.experiments{ei}
    
    mkdir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir ]);
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll.mat',expdir,exp,timeUnit);
    load(f,'gdata','keptvox');
    keptvox_wholeBrain_L=keptvox;
    gdata_listener=gdata;
    clear gdata
    
    for ri=1%:length(rois);
        roi=rois{ri};
        fr = sprintf('%s/roi_mask/%s/mat/%s_%s',expdir,froidir,roi,exp);
        load(fr,'roimask');
        
        disp(exp)
        disp(sum(roimask(keptvox_wholeBrain_L)))
        if sum(roimask(keptvox_wholeBrain_L))>0;
            clear gdata
            gdata=gdata_listener(logical(roimask(keptvox_wholeBrain_L)),:,:);
            
            if ~isempty(gdata);
                keptvox=keptvox_wholeBrain_L(ismember(keptvox_wholeBrain_L,find(roimask)));
                save([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' roi ],'gdata','keptvox');
            end
        
        end
    end
end
