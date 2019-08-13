
clc;clear;

T=20;N=5;
x0=zeros(1,2*N-1);

x0(1:N)=[3 2 5 7 3];%ʱ���theta(i)
x0(N+1:2*N-1)=[3 3 3 3];%����

lb=[0*ones(1,N) 0*ones(1,N-1)]; %�½�
ub=[T*ones(1,N) T*ones(1,N-1)];%�Ͻ�

 options = optimset('display','iter');  %������ÿ��ֵ��ʾ����
 options=optimset(options,'Algorithm','sqp');
 options=optimset(options,'tolx',1e-8);
 options = optimset(options,'GradObj','on');  %���ݶȵ���

Aeq=[ones(1,N)  zeros(1,N-1)]; beq=T;
[x,fval] = fmincon(@objgrade,x0,[],[],Aeq,beq,lb,ub,[],options)


