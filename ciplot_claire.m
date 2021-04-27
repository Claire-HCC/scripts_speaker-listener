function ciplot_claire(data,x,colour,transparency);
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% remove observation with NaN value
data=data(sum(~isnan(data),2)>0,:);

% Raymond Reynolds 24/11/06
m=nanmean(data);
ci95=ci(data,0.95);
lower=m-ci95/2;
upper=m+ci95/2;

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<3
    colour='b';
end

if nargin<2
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

h=fill([x fliplr(x)],[upper fliplr(lower)],colour)
set(h,'facealpha',transparency,'LineStyle','none');
hold on
plot(x,m,'color',colour,'linewidth',2);
hold off
