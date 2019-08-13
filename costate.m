
function dy = costate(t,y,tt,state1,state2,state3,state4,state5,taui)

%state1=S_i,state2=I_i,y=[S_i I_i]',已知y在tt=[0 1 2 3 4]时刻的取值求在t=2.5时刻的值y
E=10; R=1000; k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;a=0.01;b=0.02;
%E=10; R=1000;k1=0.5;k2=0.03;k3=3;k4=0.01;k5=0.3;k6=0.1;k7=7.2;k8=10;k9=0.9;

 
dy=zeros(5,1);
y1=interp1(tt,state1,t);
y2=interp1(tt,state2,t); 
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);%离散点tt处状态差值，t时刻的值分别为y1,y2,y3,y4,y5

dy = taui*[y(1)*(k2+k3)-y(2)*k2;%协态方程
         y(2)*(k5+k4*R/((1+a*y2)*(1+a*y2)))-y(3)*k4*R/((1+a*y2)*(1+a*y2));
         y(3)*k6*y5/((1+b*y3)*(1+b*y3))-y(4)*k6*y5/((1+b*y3)*(1+b*y3))+y(5)*k6*y5/((1+b*y3)*(1+b*y3));
         y(4)*k7;
         -1+y(3)*k6*y3-y(4)*k6*y3/(1+b*y3)+y(5)*(k9+k9*y3/(1+b*y3))];
                  
end
    
%[t,y]=ode45(@名字，[1,0],末端值，[],参数)现在可以用[1,0],前面不用加-1了
%加-1的原因（见理解图）