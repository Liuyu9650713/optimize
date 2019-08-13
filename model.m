

function model()
y0=[10 0 0 0 5];              %分别为N S SR SRM M的初始值
t0=0;
tfinal=10;
[T,Y]=ode45(@ODEmRNA,[t0,tfinal],y0);
figure(1)
plot(T,Y(:,1))
xlabel('time');
ylabel('N');
hold on

figure(2)
plot(T,Y(:,2))
xlabel('time');
ylabel('S');
hold on

figure(3)
plot(T,Y(:,3))
xlabel('time');
ylabel('SR');
hold on

figure(4)
plot(T,Y(:,4))
xlabel('time');
ylabel('SRM');
hold on

figure(5)
plot(T,Y(:,5))
xlabel('time');
ylabel('M');
hold on
end

function dydt=ODEmRNA(t,y)

E=10; R=1000; k1=0.005;k2=5e-01;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9; %根据论文给定的参数的初始值
dydt=[k1*E-(k2+k3)*y(1);
      k2*y(1)-k5*y(2)-k4*y(2)*R;
      k4*y(2)*R-k6*y(5)*y(3);
      k6*y(5)*y(3)-k7*y(4);
      k8-k9*y(5)-k6*y(5)*y(3)];
end



