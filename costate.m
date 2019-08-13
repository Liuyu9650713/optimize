
function dy = costate(t,y,tt,state1,state2,state3,state4,state5,taui)

%state1=S_i,state2=I_i,y=[S_i I_i]',��֪y��tt=[0 1 2 3 4]ʱ�̵�ȡֵ����t=2.5ʱ�̵�ֵy
E=10; R=1000; k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;a=0.01;b=0.02;
%E=10; R=1000;k1=0.5;k2=0.03;k3=3;k4=0.01;k5=0.3;k6=0.1;k7=7.2;k8=10;k9=0.9;

 
dy=zeros(5,1);
y1=interp1(tt,state1,t);
y2=interp1(tt,state2,t); 
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);%��ɢ��tt��״̬��ֵ��tʱ�̵�ֵ�ֱ�Ϊy1,y2,y3,y4,y5

dy = taui*[y(1)*(k2+k3)-y(2)*k2;%Э̬����
         y(2)*(k5+k4*R/((1+a*y2)*(1+a*y2)))-y(3)*k4*R/((1+a*y2)*(1+a*y2));
         y(3)*k6*y5/((1+b*y3)*(1+b*y3))-y(4)*k6*y5/((1+b*y3)*(1+b*y3))+y(5)*k6*y5/((1+b*y3)*(1+b*y3));
         y(4)*k7;
         -1+y(3)*k6*y3-y(4)*k6*y3/(1+b*y3)+y(5)*(k9+k9*y3/(1+b*y3))];
                  
end
    
%[t,y]=ode45(@���֣�[1,0],ĩ��ֵ��[],����)���ڿ�����[1,0],ǰ�治�ü�-1��
%��-1��ԭ�򣨼����ͼ��