%% Task 3
load('task_3_temp_check.mat','T');
T_room = mean(T(1,1:2853))
dT_room = std(T(1,1:2853))
clear T;

%% Task 4
load('task_4_temp_check','T');
T_sat_1 = mean(T(1,8500:11405))
dT_sat_1 = std(T(1,8500:11405))
clear T;

load('task_4_1_temp_check','T');
T_sat_2 = mean(T(1,6000:12000))
dT_sat_2 = std(T(1,6000:12000))
clear T;

%% Task 5-6 (calcolo H0, T1, T2)
load('task_5.mat');
t = 0.6*(1:length(T));
Vpulse = 3.5;
Vstep = 3;
eps = 18;
T = T - T_room * ones(1, length(T));
H0_1 = (T_sat_1-T_room)/Vstep
H0_2 = (T_sat_2-T_room)/Vstep 
dH0_1 = sqrt(dT_sat_1^2 + dT_room^2)/Vstep 
dH0_2 = sqrt(dT_sat_2^2 + dT_room^2)/Vstep 


%% Trovare il punto di massimo
figure();
pp = plot(t,T,'linewidth',2);

%% coordinate del punto di massimo
rp = 2.62816; % trovare a occhio e inserire a mano
tp = 27.6;

% Funzione rp = f(T12)
% vettore con T12
T12 = [1,0:0.1:400]; 
% funzione rp = f(T12)
rp_func = Vpulse*H0_1.*exp(-tp./T12).*(exp(eps./T12) - 1);
% plot della funzione con la retta di intersezione
figure();
plot(T12,rp_func,'linewidth',2);
hold on
plot(T12,rp.*ones(1, length(T12)),'linewidth',2);

% Trovare a occhio T1 e T2
T1 = 3.941;
T2 = 188.304;

%% Task 7
P=(1/3)*((T1+T2)^2)/(T1*T2)-1
I=(1/27)*((T1+T2)^3)/((T1*T2)^2)
Kp=P/H0_1
Ki=I/H0_1