g=load('DataGPS.txt');
trj=zeros(size(g,1),4);
trj(:,1)=g(:,1) - g(1,1);
trj(:,2)=fix(g(:,2)/100);
trj(:,2)=trj(:,2) + sign(g(:,2)).*(g(:,2)/100-trj(:,2))*100/60;

trj(:,3)=fix(g(:,4)/100);
trj(:,3)=trj(:,3) + sign(g(:,4)).*(g(:,4)/100-trj(:,3))*100/60;

trj(:,4)= g(:,9)+g(:,11);