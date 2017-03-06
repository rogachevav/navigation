a=load('RAV_EO');
[Pxx,F]=pwelch(a(:,1),1024,400,1024,50);
plot(F,Pxx);
set(gca,'XLim',[0 12],'YLim',[0 0.01]);