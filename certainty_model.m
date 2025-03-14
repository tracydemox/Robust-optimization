clc
clear 

%% Constant and Parameter matrix setting
f = [400, 414, 326]';
a = [18, 25, 20]';
C = [22, 33, 24;
    33, 23, 30;
    20, 25, 27];
b=[22, 33, 24, 33, 23, 30, 20, 25, 27]';
d_ = [206, 274, 220]';
% du = [40, 40, 40]';
gamma = [1.8, 1.2]';



%% Variable statement
MaxIter=10; % Max iteration
x=sdpvar(3,3,'full');
y=binvar(3,1);
z=sdpvar(3,1);
d=sdpvar(3,1);
g=sdpvar(3,1);


%% 目标函数
obj = f'*y + a'*z + sum(C.*x,'all');

C =[];
C = [C, z>=0, x>=0, 0<=g<=1];
C = [C, z<=800*y];
C = [C, sum(x,1)<=z', sum(x,2)>=d];
C = [C, d == d_ + 40*g];
ops=sdpsettings('solver','gurobi','verbose',0);
sol = optimize(C, obj, ops);

%% 结果
obj=value(obj)
y=value(y)
x=value(x)
z=value(z)
g=value(g)
d=value(d)