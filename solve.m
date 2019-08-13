
clc;clear;

T=20;N=5;
x0=zeros(1,2*N-1);

x0(1:N)=[3 2 5 7 3];%时间段theta(i)
x0(N+1:2*N-1)=[3 3 3 3];%控制

lb=[0*ones(1,N) 0*ones(1,N-1)]; %下界
ub=[T*ones(1,N) T*ones(1,N-1)];%上界

 options = optimset('display','iter');  %迭代的每个值显示出来
 options=optimset(options,'Algorithm','sqp');
 options=optimset(options,'tolx',1e-8);
 options = optimset(options,'GradObj','on');  %用梯度迭代

Aeq=[ones(1,N)  zeros(1,N-1)]; beq=T;
[x,fval] = fmincon(@objgrade,x0,[],[],Aeq,beq,lb,ub,[],options)


