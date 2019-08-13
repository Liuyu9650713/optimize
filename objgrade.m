
function [obj_f1,grad_G]=objgrade(x)
E=10; R=1000; k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;a=0.01;b=0.02;
%E=10; R=1000;k1=0.5;k2=0.03;k3=3;k4=0.01;k5=0.3;k6=0.1;k7=7.2;k8=10;k9=0.9;
P=90;

T=20;N=5;
tau=zeros(1,N);
C=zeros(1,N-1);

tau(1:N)=x(1:N);
C(1:N-1)=x(N+1:2*N-1);

y0 = [0 0 0 0 5];%状态的初始值
ya=zeros(N,5);yb=zeros(N,5);%把每一段初始值和结束时刻的状态值记录下来 ya:初始  yb:结束
y=[];t=[];%空矩阵开始默认是零矩阵，对应求解出来的时间点，状态

tn1=zeros(1,N); %每一段上离散点的个数取长之后构成的行向量，如tn1=[3,5,3,...]第一段上3个离散点，第2段上5个离散点，第3段上3个离散点...
tn2=zeros(1,N+1);%2*n+2个时间点上对应的离散点数构成的行向量,如tn2=[0,3,8,11,...]从第2个点开始等于tn1的累加求和
int_L=[]; 
 for i=1:N;%第一部分：调用状态函数，求解每一点的状态
        
         tspan=[i-1:0.1:i];
         [tt,yy] = ode45(@state,tspan,y0,[],tau(i));%tt是第i段上离散时间点的列向量，yy状态矩阵，第一列是S的数量，第二列是I的数量，形式yy(?,:)，其中？表示第几行，对应离散时刻，：表示第几个状态
         tn1(i)=length(tt);%每一段变换到[0 1]上，第i段离散时间点取长即点的数目
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:); %（行向量）把状态在0和1时刻的值记录下来,':'代表第几个状态,'i'代表第i个时间段，1代表0时刻，end代表1时刻，最后累计变成矩阵，第一列S，第二列I
         t=[t;tt]; %每一次把离散时间点的列向量累计起来，构成列向量t
         y=[y; yy]; %  y is a length(tt)*2 vector; 把每次的状态累计起来,第一列S，第二列I
          
         yya=yy(:,5);  % 目标函数原先积分号里面的部分
         int_LL=tau(i)*trapz(tt,yya);
         
         int_L=[int_L int_LL];
         
         if i<N
             y0=[C(i)*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)];
         end
 end

  tn2(2:end)=cumsum(tn1); %第2个时间点t1到最后一个时间点tN等于每一段的离散点的个数累加求和,如tn2=[0,3,8,11,...，20,23,25]从第2个点开始等于tn1的累加求和
 
  obj_f1=sum(P*E*C)+sum(int_L);
  
  grad_G=zeros(1,2*N-1);%定义梯度规则
  grad_tau=zeros(1,N);
  grad_C=zeros(1,N-1);   %定义各个参数梯度定义规则大小
  
 %下面开始调用协态函数求解协态
 
 
lam_t=[]; lamda_y=[];%空矩阵开始默认是零矩阵，对应求解出来的时间点，协态
lama=zeros(N,5);
lamb=zeros(N,5);
tn3=zeros(1,N); %每一段上离散点的个数取长之后构成的行向量，同tn1
tn4=zeros(1,N+1);%2*n+2个时间点上对应的离散点数构成的行向量,同tn2

inv_lamdai0=[0 0 0 0 0];%协态的在最后一段时间点1上的初始值

 for i=1:N;
     
     %i对应第N-i+1个时间段先求第5段，然后第4段，...
     
     ti=t(tn2(N-i+1)+1:tn2(N-i+2));%每一段上的时间点构成的列向量,（对应t里的元素），如果最后一段是两个时间点0.2,0.9，则t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',倒求
     yi=y(tn2(N-i+1)+1:tn2(N-i+2),:);%把每一段上每个点的状态记录下来 （对应y里的元素），‘：’代表状态1,状态2,第一列代表S，第2列代表I，倒求
     
     inv_tspan=[N-i+1:-0.1:N-i];

     [inv_tti,inv_lamdai]=ode45(@(t,y)costate(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5),tau(N-i+1)),inv_tspan, inv_lamdai0);
     %inv_tti是第N-i+1段上离散时间点列向量，inv_lamdai是协态的矩阵，第一列是S的协态，第二列是I的协态，形式inv_lamdai(?,:)，其中？表示第几个时刻，：表示第几个协态

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai'); %时间，协态翻转
 
     lama(N-i+1,:)=lamdai(:,1)'; lamb(N-i+1,:)=lamdai(:,end)';
     tn3(N-i+1)=length(tti);%每一段上对应协态离散化的时间点取长,tn3倒取
     
     lam_t=[tti lam_t];   %时间计加到一块
     lamda_y=[lamdai lamda_y];   %协态累计到一块
              
         if i<N  
          
         inv_lamdai0=[lama(N-i+1,1) lama(N-i+1,2) lama(N-i+1,3) lama(N-i+1,4) lama(N-i+1,5)];% 协态跳跃
               
         end
 end

 
  tn4(2:end)=cumsum(tn3); %同tn2
 
 %下面把状态和协态进行标准化
 
 for i=1:N;
     
dt=0.02;  
stan_t=i-1:dt:i;  %把协态和状态值统一在相同的离散点上

%下面先把状态标准化
stan_y1=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),1),stan_t);
stan_y2=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),2),stan_t);
stan_y3=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),3),stan_t);
stan_y4=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),4),stan_t);
stan_y5=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),5),stan_t);

%把协态标准化
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
 




