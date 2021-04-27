clear all
set_parameters;
% framewise displacement (FD): summing the absolute values of the differentials of the six parameters.
for ei=11:12;%:4;
    exp=experiments{ei};
    displacements=nan([48 1]);
    movements=nan([48 1 6]);
    
    for si=1:48;
        f=sprintf('%s/%s/movements/listeners/sub-%02d_task-%s_bold_confounds.tsv',expdir,exp,si,strrep(exp,'_fmriprep_noConfoundRegression',''));
        if exist(f)~=0;
        f=tdfread(f);
        displacement=str2num(f.FramewiseDisplacement(2:end,:));
        displacements(si,1:length(displacement))=displacement;
        movement=[f.X f.Y f.Z f.RotX f.RotY f.RotZ];
        movements(si,1:length(movement),:)=movement;
        end
    end
    outf=sprintf('%s/%s/movements/listeners/displacements.mat',expdir,exp);
    save(outf,'displacements','movements');
    find(max(displacements,[],2)>3)
    find(max(max(movements,[],2),[],3)>3)
end
% piement 22, translation along z-axis > 3mm

for ei=11:12%:2;%:4;
    exp=experiments{ei};
    f=dir([expdir exp '/movements/listeners/*tsv'])
    f=tdfread([f.folder '\' f.name]);
    displacement=str2num(f.FramewiseDisplacement(2:end,:));
    movement=[f.X f.Y f.Z f.RotX f.RotY f.RotZ];
    figure;
    subplot(2,1,1);
    plot(displacement);
    title({exp,'framewise-displacement'});
    grid on
    subplot(2,1,2);
    plot(movement);
    title('movements');
    grid
end






