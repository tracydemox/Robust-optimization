%% CCG algorithm for Two-stage Robust Optimization
% [1] Zeng, Bo, and Long Zhao. "Solving two-stage robust optimization problems using a column-and-constraint generation method." Operations Research Letters 41, no. 5 (2013): 457-461.
clear 
clc
warning off

%% Constant and Parameter matrix setting
f = [400, 414, 326]';
a = [18, 25, 20]';
C = [22, 33, 24;
    33, 23, 30;
    20, 25, 27];
b=[22, 33, 24, 33, 23, 30, 20, 25, 27]';
dl = [206, 274, 220]';
du = [40, 40, 40]';
G1=[1,1,1,0,0,0,0,0,0;
    0,0,0,1,1,1,0,0,0;
    0,0,0,0,0,0,1,1,1];
G2=[1,0,0,1,0,0,1,0,0;
    0,1,0,0,1,0,0,1,0;
    0,0,1,0,0,1,0,0,1];
BigM=1e4;

%% Variable statement
MaxIter=10; % Max iteration
% 第一阶段决策变量
% x=sdpvar(9,MaxIter);
x=sdpvar(3,3,MaxIter,'full');
% 第二阶段决策变量
y=binvar(3,1);
z=sdpvar(3,1);
% 不确定量
d=sdpvar(3,1);
g=sdpvar(3,1);
% 拉格朗日乘子
lamda=sdpvar(3,1);
pai=sdpvar(3,1);
% miu=sdpvar(9,1);
miu=sdpvar(3,3,'full');
% 辅助变量
% v=binvar(9,1);
v=binvar(3,3,'full');
s=binvar(3,1);
q=binvar(3,1);
eta=sdpvar(1);

%% CCG algorithm
LB=-inf;
UB=inf;
Epsilon=0.01;

Obj_MP = f'*y + a'*z + eta;
Cons_MP =[z<=800*y, z>=0, eta>=sum(C.*x(:,:,1),'all'),  sum(z)>=772 , x(:,:,1)>=0 ];
% Cons_MP =[z<=800*y, z>=0, eta>=b'*x(:,1),  sum(z)>=772,  x(:,1)>=0 ];

ops=sdpsettings('solver','gurobi','verbose',0);
% ops.gurobi.MIPGap = 1e-2; 
% ops.gurobi.TimeLimit = 5;      % 设置求解时间为10s
Cons_SP = [];
tic;

for k=1:5
    % Solve MP
    sol_MP=optimize(Cons_MP,Obj_MP,ops);
    s_y=value(y);
    s_z=value(z);
    s_eta=value(eta);
    LB=value(Obj_MP);

    % Solve SP with kkt
%     Obj_SP = b'*x(:,k);
    Obj_SP = sum(C.*x(:,:,k),'all');
    Cons_SP = [Cons_SP, pai>=0, lamda>=0, miu>=0];
%     Cons_SP = [Cons_SP, G1*x(:,k)<=s_z, G2*x(:,k) >= dl + 40*g, x(:,k)>=0];
Cons_SP = [Cons_SP, sum(x(:,:,k),2)<=s_z, sum(x(:,:,k))' >= dl + 40*g, x(:,k)>=0];

    % 对偶可行性条件
%     Cons_SP = [Cons_SP,  b + G1'*lamda - G2'*pai - miu ==0];
for i=1:3
    for j=1:3
        Cons_SP = [Cons_SP,  C(i,j) + lamda(i) - pai(j) - miu(i,j) ==0];
        Cons_SP = [Cons_SP, miu(i,j)<=BigM*v(i,j), -x(i,j,k)>=BigM*(v(i,j)-1)];
        Cons_SP = [Cons_SP, pai<=BigM*s, dl + 40*g - sum(x(:,:,k))'>=BigM*(s-1)];
        Cons_SP = [Cons_SP, lamda<=BigM*q,  sum(x(:,:,k),2)-s_z >=BigM*(q-1)];
    end
end
    % 互补松弛条件 --- big-M法线性化
%     Cons_SP = [Cons_SP, miu<=BigM*v, -x(:,k)>=BigM*(v-1)];
%     Cons_SP = [Cons_SP, pai<=BigM*s, dl + 40*g - G2*x(:,k)>=BigM*(s-1)];
%     Cons_SP = [Cons_SP, lamda<=BigM*q,  G1*x(:,k)-s_z >=BigM*(q-1)];

    Cons_SP = [Cons_SP, sum(g)<=1.8, g(1)+g(2)<=1.2, 0<=g<=1];

    sol_SP = optimize(Cons_SP,-Obj_SP,ops);
    UB=min(UB, f'*s_y+a'*s_z+value(Obj_SP));
                                                                                                                                                                                                                             
% if sol_SP.problem == 0
%     disp('succcessful solved');
% else
%     yalmiperror(sol_SP.problem)
% end

disp(['第',num2str(k),'次迭代']);
disp(['UB: ',num2str(UB),' LB: ',num2str(LB)]);

if UB-LB<=Epsilon
    break;
end
%   Cons_MP = [Cons_MP, eta>=b'*x(:,k+1), G1*x(:,k+1)<=z, G2*x(:,k+1)>= dl + 40*value(g),  x(:,k+1)>=0];
  Cons_MP = [Cons_MP, eta>=sum(C.*x(:,:,k+1),'all'), sum(x(:,:,k+1),2)<=z, sum(x(:,:,k+1))'>= dl + 40*value(g),  x(:,:,k+1)>=0];
end

toc;