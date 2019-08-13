
function [obj_f1,grad_G]=objgrade(x)
E=10; R=1000; k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;a=0.01;b=0.02;
%E=10; R=1000;k1=0.5;k2=0.03;k3=3;k4=0.01;k5=0.3;k6=0.1;k7=7.2;k8=10;k9=0.9;
P=90;

T=20;N=5;
tau=zeros(1,N);
C=zeros(1,N-1);

tau(1:N)=x(1:N);
C(1:N-1)=x(N+1:2*N-1);

y0 = [0 0 0 0 5];%״̬�ĳ�ʼֵ
ya=zeros(N,5);yb=zeros(N,5);%��ÿһ�γ�ʼֵ�ͽ���ʱ�̵�״ֵ̬��¼���� ya:��ʼ  yb:����
y=[];t=[];%�վ���ʼĬ��������󣬶�Ӧ��������ʱ��㣬״̬

tn1=zeros(1,N); %ÿһ������ɢ��ĸ���ȡ��֮�󹹳ɵ�����������tn1=[3,5,3,...]��һ����3����ɢ�㣬��2����5����ɢ�㣬��3����3����ɢ��...
tn2=zeros(1,N+1);%2*n+2��ʱ����϶�Ӧ����ɢ�������ɵ�������,��tn2=[0,3,8,11,...]�ӵ�2���㿪ʼ����tn1���ۼ����
int_L=[]; 
 for i=1:N;%��һ���֣�����״̬���������ÿһ���״̬
        
         tspan=[i-1:0.1:i];
         [tt,yy] = ode45(@state,tspan,y0,[],tau(i));%tt�ǵ�i������ɢʱ������������yy״̬���󣬵�һ����S���������ڶ�����I����������ʽyy(?,:)�����У���ʾ�ڼ��У���Ӧ��ɢʱ�̣�����ʾ�ڼ���״̬
         tn1(i)=length(tt);%ÿһ�α任��[0 1]�ϣ���i����ɢʱ���ȡ���������Ŀ
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:); %������������״̬��0��1ʱ�̵�ֵ��¼����,':'����ڼ���״̬,'i'�����i��ʱ��Σ�1����0ʱ�̣�end����1ʱ�̣�����ۼƱ�ɾ��󣬵�һ��S���ڶ���I
         t=[t;tt]; %ÿһ�ΰ���ɢʱ�����������ۼ�����������������t
         y=[y; yy]; %  y is a length(tt)*2 vector; ��ÿ�ε�״̬�ۼ�����,��һ��S���ڶ���I
          
         yya=yy(:,5);  % Ŀ�꺯��ԭ�Ȼ��ֺ�����Ĳ���
         int_LL=tau(i)*trapz(tt,yya);
         
         int_L=[int_L int_LL];
         
         if i<N
             y0=[C(i)*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)];
         end
 end

  tn2(2:end)=cumsum(tn1); %��2��ʱ���t1�����һ��ʱ���tN����ÿһ�ε���ɢ��ĸ����ۼ����,��tn2=[0,3,8,11,...��20,23,25]�ӵ�2���㿪ʼ����tn1���ۼ����
 
  obj_f1=sum(P*E*C)+sum(int_L);
  
  grad_G=zeros(1,2*N-1);%�����ݶȹ���
  grad_tau=zeros(1,N);
  grad_C=zeros(1,N-1);   %������������ݶȶ�������С
  
 %���濪ʼ����Э̬�������Э̬
 
 
lam_t=[]; lamda_y=[];%�վ���ʼĬ��������󣬶�Ӧ��������ʱ��㣬Э̬
lama=zeros(N,5);
lamb=zeros(N,5);
tn3=zeros(1,N); %ÿһ������ɢ��ĸ���ȡ��֮�󹹳ɵ���������ͬtn1
tn4=zeros(1,N+1);%2*n+2��ʱ����϶�Ӧ����ɢ�������ɵ�������,ͬtn2

inv_lamdai0=[0 0 0 0 0];%Э̬�������һ��ʱ���1�ϵĳ�ʼֵ

 for i=1:N;
     
     %i��Ӧ��N-i+1��ʱ��������5�Σ�Ȼ���4�Σ�...
     
     ti=t(tn2(N-i+1)+1:tn2(N-i+2));%ÿһ���ϵ�ʱ��㹹�ɵ�������,����Ӧt���Ԫ�أ���������һ��������ʱ���0.2,0.9����t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',����
     yi=y(tn2(N-i+1)+1:tn2(N-i+2),:);%��ÿһ����ÿ�����״̬��¼���� ����Ӧy���Ԫ�أ�������������״̬1,״̬2,��һ�д���S����2�д���I������
     
     inv_tspan=[N-i+1:-0.1:N-i];

     [inv_tti,inv_lamdai]=ode45(@(t,y)costate(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5),tau(N-i+1)),inv_tspan, inv_lamdai0);
     %inv_tti�ǵ�N-i+1������ɢʱ�����������inv_lamdai��Э̬�ľ��󣬵�һ����S��Э̬���ڶ�����I��Э̬����ʽinv_lamdai(?,:)�����У���ʾ�ڼ���ʱ�̣�����ʾ�ڼ���Э̬

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai'); %ʱ�䣬Э̬��ת
 
     lama(N-i+1,:)=lamdai(:,1)'; lamb(N-i+1,:)=lamdai(:,end)';
     tn3(N-i+1)=length(tti);%ÿһ���϶�ӦЭ̬��ɢ����ʱ���ȡ��,tn3��ȡ
     
     lam_t=[tti lam_t];   %ʱ��Ƽӵ�һ��
     lamda_y=[lamdai lamda_y];   %Э̬�ۼƵ�һ��
              
         if i<N  
          
         inv_lamdai0=[lama(N-i+1,1) lama(N-i+1,2) lama(N-i+1,3) lama(N-i+1,4) lama(N-i+1,5)];% Э̬��Ծ
               
         end
 end

 
  tn4(2:end)=cumsum(tn3); %ͬtn2
 
 %�����״̬��Э̬���б�׼��
 
 for i=1:N;
     
dt=0.02;  
stan_t=i-1:dt:i;  %��Э̬��״ֵ̬ͳһ����ͬ����ɢ����

%�����Ȱ�״̬��׼��
stan_y1=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),1),stan_t);
stan_y2=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),2),stan_t);
stan_y3=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),3),stan_t);
stan_y4=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),4),stan_t);
stan_y5=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),5),stan_t);

%��Э̬��׼��
stan_lam1=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(1,tn4(i)+1:tn4(i+1)),stan_t);
stan_lam2=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(2,tn4(i)+1:tn4(i+1)),stan_t);
stan_lam3=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(3,tn4(i)+1:tn4(i+1)),stan_t);
stan_lam4=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(4,tn4(i)+1:tn4(i+1)),stan_t);
stan_lam5=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(5,tn4(i)+1:tn4(i+1)),stan_t);

stan_ylam=stan_y5+stan_lam1.*(k1*E-(k2+k3)*stan_y1)+stan_lam2.*(k2*stan_y1-k5*stan_y2-k4*stan_y2*R/(1+a*stan_y2))+stan_lam3.*(k4*R*stan_y2/(1+a*stan_y2)-k6*stan_y3.*stan_y5/(1+b*stan_y3))+...
          stan_lam4.*(k6*stan_y3.*stan_y5/(1+b*stan_y3)-k7*stan_y4)+stan_lam5.*(k8-k9*stan_y5-k6*stan_y3.*stan_y5/(1+b*stan_y3));
     
          
             if i<N
           
             grad_C(i)=P*E+E*lama(i+1,1);
             
             end 
             
   grad_tau(i)=trapz(stan_t,stan_ylam);         
 
 end
 
 grad_G=[grad_tau grad_C];
 
 
end
 




