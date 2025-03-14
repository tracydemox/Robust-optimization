%% CCG algorithm for Two-stage Robust Optimization
% [1] Zeng, Bo, and Long Zhao. "Solving two-stage robust optimization problems using a column-and-constraint generation method." Operations Research Letters 41, no. 5 (2013): 457-461.
clear 
clc
%% Constant and Parameter matrix setting
f = [400, 414, 326]';
a = [18, 25, 20]';
C = [22, 33, 24;
    33, 23, 30;
    20, 25, 27];
b=[22, 33, 24, 33, 23, 30, 20, 25, 27]';
dl = [206, 274, 220]';
du = [40, 40, 40]';
G=[-1,-1,-1,0,0,0,0,0,0;
  0,0,0,-1,-1,-1,0,0,0;
  0,0,0,0,0,0,-1,-1,-1;
  1,0,0,1,0,0,1,0,0;
  0,1,0,0,1,0,0,1,0;
  0,0,1,0,0,1,0,0,1];
G1=[-1,-1,-1,0,0,0,0,0,0;
    0,0,0,-1,-1,-1,0,0,0;
    0,0,0,0,0,0,-1,-1,-1];
G2=[1,0,0,1,0,0,1,0,0;
    0,1,0,0,1,0,0,1,0;
    0,0,1,0,0,1,0,0,1];
BigM=1e5;

%% Variable statement
MaxIter=10; % Max iteration
x=sdpvar(9,MaxIter);
y=binvar(3,1);
z=sdpvar(3,1);
% 不确定量
d=sdpvar(3,1);
g=sdpvar(3,1);
% g=binvar(3,1);
% 对偶变量
alpha=sdpvar(3,1);
beta=sdpvar(3,1);
% 辅助变量
w=sdpvar(3,1);
eta=sdpvar(1);

%% CCG algorithm
LB=-inf;
UB=inf;
Epsilon=0.01;

Obj_MP2 = f'*y + a'*z + eta;
Cons_MP2 =[z<=800*y, z>=0, eta>=b'*x(:,1), sum(z)>=772, x(:,1)>=0 ];

ops=sdpsettings('solver','gurobi','verbose',0);
ops.gurobi.MIPGap = 1e-2; 
ops.gurobi.TimeLimit = 5;      % 设置求解时间为10s
Cons_SP2 = [];


for k=1:5
  % Solve MP2

    sol_MP2=optimize(Cons_MP2,Obj_MP2,ops);
    s_y=value(y);
    s_z=value(z);
    s_eta=value(eta);
    LB=value(Obj_MP2);

 % Solve SP2
    Obj_SP2 = - s_z'*alpha + dl'*beta + 40*g'*beta ;
%     Obj_SP2 = -s_z'*alpha + dl'*beta + du'*w ;
    Cons_SP2 = [Cons_SP2, alpha>=0, beta>=0 ];
%     for i=1:3
%         for j=1:3
%             Cons_SP2 = [Cons_SP2, -alpha(i) + beta(j) <= C(i,j) ];
%         end
%     end
     Cons_SP2 = [Cons_SP2, G1'*alpha + G2'*beta <= b ];
%     Cons_SP2 = [Cons_SP2, beta-BigM*(1-g)<=w<=beta, 0<=w<=BigM*g];
%     Cons_SP2 = [Cons_SP2, sum(g)<=2];

 Cons_SP2 = [Cons_SP2, sum(g)<=1.8, g(1)+g(2)<=1.2, 0<=g<=1];

    sol_SP2 = optimize(Cons_SP2,-Obj_SP2,ops);
    UB=min(UB, f'*s_y+a'*s_z+value(Obj_SP2));

% if sol_SP2.problem == 0
%     disp('succcessful solved');
% else
%     disp('error');
%     yalmiperror(sol.problem)
% end
disp(['第',num2str(k),'次迭代']);
disp(['UB: ',num2str(UB),' LB: ',num2str(LB)])

if UB-LB<=Epsilon
    break;
end

Cons_MP2 = [Cons_MP2, eta>=b'*x(:,k+1), G1*x(:,k+1)>=-z, G2*x(:,k+1)>= dl+40*value(g),  x(:,k+1)>=0];
end
 