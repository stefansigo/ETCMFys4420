clear all;
clc;

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
load('task_5_temp_check_1.mat');
Vpulse = 3.5;
Vstep = 3;
eps = 18;
T = T - T_room * ones(1, length(T));
H0_1 = (T_sat_1-T_room)/Vstep
H0_2 = (T_sat_2-T_room)/Vstep 
dH0_1 = sqrt(dT_sat_1^2 + dT_room^2)/Vstep
dH0_2 = sqrt(dT_sat_2^2 + dT_room^2)/Vstep

% Trovare massimo
figure();
plot(cm(1,1:300),T, '.', 'markersize', 10); %,'color', '[0.25 0.9 0.25]');
grid on;

%% coordinate del punto di massimo
rp = 2.62816 % trovare a occhio e inserire a mano
tp = 27.05

% Funzione rp = f(T12)
% vettore con T12
T12 = [1,0:0.1:400]; 
% funzione rp = f(T12)
rp_func = Vpulse*H0_1.*exp(-tp./T12).*(exp(eps./T12) - 1);
% plot della funzione con la retta di intersezione
figure();
plot(T12,rp_func, '.', 'markersize', 10);
hold on
plot(T12,rp.*ones(1, length(T12)), '.', 'markersize', 10);
grid on;

%% Trovare a occhio T1 e T2
T1 = 190.81
T2 = 3.6977

%% Task 7
P=(1/3)*((T1+T2)^2)/(T1*T2)-1
I=(1/27)*((T1+T2)^3)/((T1*T2)^2)
Kp=P/H0_1
Ki=I/H0_1

% dP = (1/3)*sqrt((0.0001)^2 + (0.0001))*P;
dP = 0.01;
dKp = sqrt(dP^2 + dH0_1^2)*Kp
% dI = (1/27)*sqrt((0.0001)^3 + (0.0001)^2)*I;
dI = 0.0002;
dKi = sqrt(dI^2 + dH0_1^2)*Ki

%% Task 12
delete T cm Tref
load('task_12_pid_control.mat');

% Rapporto sperimentale prima reference temperature
Tas1 = mean(T(327:483)) - T_room;
dTas1 = sqrt(std(T(327:483))^2+dT_room^2);
Tref1 = Tref(327) - T_room
dTref1 = dT_room;
Rapp1 = Tas1/Tref1
dRapp1 = sqrt((dTas1/Tref1)^2 + ((Tas1/Tref1^2)*dTref1)^2)

% Rapporto sperimentale seconda reference temperature
Tas2 = mean(T(1183:1500)) - T_room;
dTas2 = sqrt(std(T(1183:1500))^2+dT_room^2);
Tref2 = Tref(1400) - T_room;
dTref2 = dT_room;
Rapp2 = Tas2/Tref2
dRapp2 = sqrt((dTas2/Tref2)^2 + ((Tas2/Tref2^2)*dTref2)^2)

% Rapporto teorico
RapTeo = Kp/(H0_1^(-1) + Kp)
dRapTeo = sqrt((H0_1/(1+Kp*H0_1)^2*dKp)^2+(Kp/(1+Kp*H0_1)^2*dH0_1)^2)

% Residui per compatibilità
Res1 = abs(Rapp1 - RapTeo);
dRes1 = sqrt(dRapp1^2+dRapTeo^2);
Res2 = abs(Rapp2 - RapTeo);
dRes2 = sqrt(dRapp2^2+dRapTeo^2);
