r=load('RR1.txt');
t=r(:,1);
k=180/pi;
kd=250;
psi=r(:,4); gamma=r(:,3); tetta=r(:,2);
ind=load('index.txt');
kq=10*ones([1,size(t)]);
kj=kq';
figure(1);
grid on;
hold on;
plot (t,unwrap(gamma).*k,'g','LineWidth',2);
plot(t,unwrap(psi).*k,'r','LineWidth',2);

plot (t,tetta.*k ,'b','LineWidth',2);
% plot(ind(:,1),ind(:,2).*10,'*');
% plot(ind(:,1),ind(:,3).*5,'*');
% plot(ind(:,1),ind(:,4),'*');
% plot(ind(:,1),ind(:,5).*3,'*');
% plot(ind(:,1),ind(:,6).*15,'*');
hold off;

figure(2);
grid on;
hold on;
plot(t,unwrap(r(:,5)),'r','LineWidth',2)
plot(t,unwrap(r(:,6)),'g','LineWidth',2)
plot(t,unwrap(r(:,7)),'b','LineWidth',2)
hold off;

figure(3);
grid on;
hold on;
plot(t,r(:,8).*k,'r','LineWidth',2)
plot(t,r(:,9).*k,'g','LineWidth',2)
hold off;
r(end,8)*k 
r(end,9)*k