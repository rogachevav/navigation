f=load('Filter.txt');
t=f(:,1);

figure(1);
grid on;
hold on;
plot(t,f(:,2),'r','LineWidth',2)
plot(t,f(:,3),'g','LineWidth',2)
plot(t,f(:,4),'b','LineWidth',2)
hold off;

figure(4);
grid on;
hold on;
plot(t,f(:,5),'r','LineWidth',2)
plot(t,f(:,6),'g','LineWidth',2)
plot(t,f(:,7),'b','LineWidth',2)
hold off;

figure(3);
grid on;
hold on;
plot(t,f(:,8),'r','LineWidth',2)
plot(t,f(:,9),'g','LineWidth',2)
plot(t,f(:,10),'b','LineWidth',2)
hold off;