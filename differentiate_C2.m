clc;clear;

E=10; R=1000;k1=0.005;k2=5e-04;k3=3;k4=0.001;k5=0.03;k6=0.1;k7=7.2;k8=10;k9=0.9;P=9;

N=5;
tau=[4 4 4 4 4];%ʱ���theta(i)
dt=0.1;C2=0:dt:5;
C1=2;C3=2;C4=2;

ya=zeros(N,5);yb=zeros(N,5);%��ÿ�εĳ�ʼ�ͽ���ʱ�̶���¼������ya��¼��ʼ״̬��yb�������״̬

y=[];t=[];
tn1=zeros(1,N);
tn2=zeros(1,N+1);

int_L=[];
obj_f1=[];
grad_C2=[];
bbb=[];
ccc=[];


for j=1:5/dt+1;
    y0=[0 0 0 0 5];%�����ʼ״̬�ĳ�ֵ
    
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
        
                        
                        yya=yy(:,5);  % Ŀ�꺯��ԭ�Ȼ��ֺ�����Ĳ���
                        int_LL=tau(i)*trapz(tt,yya);%����ÿһ�εĻ���
                        int_L=[int_L int_LL];%��5�εĻ��ִ洢����
                      
       end
                    
        tn2(2:end)=cumsum(tn1);%tn1�������
        bb=sum(int_L);
        bbb=[bbb bb];
        obj_f=P*C2(j)+P*C1+P*C3+P*C4+sum(int_L);%�ɱ�����
        obj_f1=[obj_f1 obj_f];
                    
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
                        
         lama(N-i+1,:)=lamdai(:,1)';
         lamb(N-i+1,:)=lamdai(:,end)';
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
                        
       end
       grad_CC2=P+E*lama(3,1);
       grad_C2=[grad_C2 grad_CC2];
                    
                    
      end
      dt=0.1;
      cc=diff(obj_f1);
      ccc=[ccc cc];
      obj_fgrad=diff(obj_f1)/dt;%�����������25ά
      aaa=5/dt;
      aa=C2(1:5/dt);
                
      figure(1)
      plot(C2,obj_f1,'b','LineWidth',2)
      hold on
      figure(4)
      plot(C2,grad_C2,'ko',C2(1:5/dt),obj_fgrad,'ro','LineWidth',1)%��ɫ���ǲ��
      hold on
                
