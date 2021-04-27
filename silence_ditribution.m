clear all
% close all
set_parameters;

eis=[11 12];
edges=[0:0.1:5];
centers=edges(1:(end-1))+diff(edges)/2;
figure;
N=[];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    tb=readtable([expdir exp '/sound/' exp '_listener_silence.xlsx']);
    dur_silence=tb.dur_sec(tb.silence==1 & tb.tmax_sec<voln_listener{ei}*tr(ei));
    [N(end+1,:)  ]=histcounts( dur_silence,edges);
%    N(end,:)=N(end,:)/length(dur_silence);
   %  max( dur_silence)
    sum(dur_silence>1.5)
    %   histogram(silence,0:0.1:5);
    %  hold on;
end

figure;
bar(N', 'grouped');
set(gca,'xticklabels',centers,'fontsize',14)
legend(experiments(eis));
xlabel('duration of silent pause (sec)');
ylabel('Number of pause');
hold off
legend boxoff