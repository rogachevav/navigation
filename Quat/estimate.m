r=load('RR1.txt');
f=load('Filter.txt');
trj=zeros(146750,4);
trj(:,1)= r(100*250:end,1);

trj(:,2)=r(100*250:end,5) - f(100*250:end,5);
trj(:,3)=r(100*250:end,6) - f(100*250:end,6);
trj(:,4)=r(100*250:end,7) - f(100*250:end,7);


figure(1);
grid on;
hold on;
plot(trj(:,1),trj(:,2),'r','LineWidth',2)
plot(trj(:,1),trj(:,3),'g','LineWidth',2)
plot(trj(:,1),trj(:,4),'b','LineWidth',2)
hold off;