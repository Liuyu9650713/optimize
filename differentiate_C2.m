clc;clear;

E=10; R=1000;k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;P=9;

N=5;
tau=[4 4 4 4 4];%时间段theta(i)
dt=0.1;C2=0:dt:5;
C1=2;C3=2;C4=2;

ya=zeros(N,5);yb=zeros(N,5);%把每段的初始和结束时刻都记录下来，ya记录初始状态，yb代表结束状态

y=[];t=[];
tn1=zeros(1,N);
tn2=zeros(1,N+1);

int_L=[];
obj_f1=[];
grad_C2=[];
bbb=[];
ccc=[];


for j=1:5/dt+1;
    y0=[0 0 0 0 5];%定义初始状态的初值
    
    for i=1:N
        Tspan=[i-1:0.1:i];
        [tt,yy]=ode45(@state,Tspan,y0,[],tau(i)); %solve state, tt is the column vector, yy is the matrix(column:S,I,P,V)
        tn1(i)=length(tt); %the number of the points
        ya(i,:)=yy(1,:); yb(i,:)=yy(end,:);
        t=[t;tt];  y=[y;yy];
        if i<N
        switch i
            case 1
                y0=[C1*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)];break
                case 2
                    y0=[C2(j)*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)];break
                    
                case 3
                        y0=[C3*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)];break
                case 4
                         y0=[C4*E+yb(i,1)  yb(i,2)  yb(i,3)  yb(i,4)  yb(i,5)]; break
                            
        end
        end
        
                        
                        yya=yy(:,5);  % 目标函数原先积分号里面的部分
                        int_LL=tau(i)*trapz(tt,yya);%计算每一段的积分
                        int_L=[int_L int_LL];%把5段的积分存储起来
                      
       end
                    
        tn2(2:end)=cumsum(tn1);%tn1点数求和
        bb=sum(int_L);
        bbb=[bbb bb];
        obj_f=P*C2(j)+P*C1+P*C3+P*C4+sum(int_L);%成本，数
        obj_f1=[obj_f1 obj_f];
                    
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
                        
         lama(N-i+1,:)=lamdai(:,1)';
         lamb(N-i+1,:)=lamdai(:,end)';
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
                        
       end
       grad_CC2=P+E*lama(3,1);
       grad_C2=[grad_C2 grad_CC2];
                    
                    
      end
      dt=0.1;
      cc=diff(obj_f1);
      ccc=[ccc cc];
      obj_fgrad=diff(obj_f1)/dt;%行向量，差导，25维
      aaa=5/dt;
      aa=C2(1:5/dt);
                
      figure(1)
      plot(C2,obj_f1,'b','LineWidth',2)
      hold on
      figure(4)
      plot(C2,grad_C2,'ko',C2(1:5/dt),obj_fgrad,'ro','LineWidth',1)%红色的是差分
      hold on
                
