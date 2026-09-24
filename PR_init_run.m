%% clear
clear
close all
clc
%% initialize
Vdc = 800;
V2 = 400;
Vpeak = V2*sqrt(2);
Tsim = 0.5;
Tstep = 0.06;
Tdist = 0.3;
TQ  = 0.08;
L = 1.3e-3;
R = 1e-2;
Rcap = 50e-3;
f = 50;
w = 2*pi*f;
wc = 0.6283;
max_i = 30;

Ki = 120;
Kp = 1;

tau = 1e-3;
kipq = 1.5*V2*sqrt(2);
dr = sqrt(2)/2;
tc = tau;
tp = 10*tc;
kpp = tc/(kipq*tp);
kip = 1/(kipq*tp);
sogi_kp = 0.5;
sogi_ki = 30;
kp = 20;
ki = 1000;
eng_kp = 20;
eng_ki = 1000;

syms V1 delta
X = (1.3e-3)*2*pi*50;
C = 2e-3;
P = 10e3;
Q = 0;
Vdc = 800;
Eref1 = 3*(0.5)*C*Vdc^2;
Eref2 = 0;
Eref3 = 0;
eqP = (V1*V2*sin(delta))/X == P;
eqQ = (V1*V2*cos(delta)-V2^2)/X == Q;
res  = solve([eqP,eqQ],[V1,delta]);
Vm = max(double(res.V1));
angle = min(abs(double(res.delta)));
P1ref = 0;
Q1ref = -15e3;
P2ref = 0;
P3ref = 0;
Ptrf = [1 1 1;1 -1 0; 1 0 -1;]; 

syms ialmag ialph ibemag ibeph
% valmag = 565.69;
% valph = 1.572;
% vbemag = 565.69;
% vbeph = -3.143;
% P1 = (3/2)*(valmag*ialmag/2)*cos(valph - ialph) + (3/2)*(vbemag*ibemag/2)*cos(vbeph - ibeph) == P1ref;
% Q1 = (3/2)*(valmag*ialmag/2)*sin(valph - ialph) + (3/2)*(vbemag*ibemag/2)*sin(vbeph - ibeph) == Q1ref;
% P2 = (3/4)*(valmag*ialmag/2)*cos(valph - ialph) - (3/4)*(vbemag*ibemag/2)*cos(vbeph - ibeph) + ...
%      (sqrt(3)/4)*(ialmag*vbemag/2)*cos(vbeph - ialph) + (sqrt(3)/4)*(ibemag*valmag/2)*cos(valph - ibeph) == P2ref;
% P3 = (3/4)*(valmag*ialmag/2)*cos(valph - ialph) - (3/4)*(vbemag*ibemag/2)*cos(vbeph - ibeph) - ...
%      (sqrt(3)/4)*(ialmag*vbemag/2)*cos(vbeph - ialph) - (sqrt(3)/4)*(ibemag*valmag/2)*cos(valph - ibeph) == P3ref;
% 
% sol = solve([P1 Q1 P2 P3],[ialmag ialph ibemag ibeph]);
% ialmag_res= double(sol.ialmag)
% ialph_res = double(sol.ialph)
% ibemag_res = double(sol.ibemag)
% ibeph_res = double(sol.ibeph)
%% run sim
out = sim('ThreePhase_Energy.slx');
%% Results / Plotting
res_PQ = out.PQ.signals.values;
res_curr = out.currents.signals.values;
res_volt = out.volts.signals.values;
res_SOGI = out.SOGI.signals.values;
res_eng = out.energy.signals.values;
e_eng = out.errors.signals.values(:,3:end);

alph = res_SOGI(:,1);
almag = res_SOGI(:,2);
beph = res_SOGI(:,3);
bemag = res_SOGI(:,4);
P_res = res_PQ(:,1);
Q_res = res_PQ(:,2);
abc_power = res_PQ(:,3);
Ptot_ref = res_PQ(:,4);
Qtot_ref = res_PQ(:,5);
e_ia = out.errors.signals.values(:,1);
e_ib = out.errors.signals.values(:,2);
ref_a = res_curr(:,4);
ref_b = res_curr(:,5);
i_a = res_curr(:,1);
i_b = res_curr(:,2);
t = out.PQ.time;
vg_al = res_volt(:,1);
vg_be = res_volt(:,2);
vc_al = res_volt(:,4);
vc_be = res_volt(:,5);
vc_abc = res_volt(:,7:9);
vg_abc = res_volt(:,10:12);
close all

% f1 = figure('Position',[200 100 550 400]);
% plot(t,i_a,'LineWidth',1.5,'DisplayName','$I_{alpha}$');
% hold on
% plot(t,ref_a,'LineWidth',1.5,'DisplayName','$I_{ref}$');
% xlabel('Time (s)')
% ylabel('Current (A)')
% title('Current vs. Error')
% legend('Interpreter','latex')

f2 = figure('Position',[200 100 550 400]);
plot(t,e_ia,'LineWidth',1.5,'DisplayName','$e_\alpha$');
hold on
plot(t,e_ib,'LineWidth',1.5,'DisplayName','$e_\beta$');
xlabel('Time (s)')
ylabel('Current (A)')
title('Current Errors')
legend('Interpreter','latex')

% f3 = figure('Position',[200 100 550 400]);
% plot(t,P_res,'LineWidth',1.5,'DisplayName','P');
% hold on
% plot(t,P1ref*ones(1,length(t)),'LineWidth',1.5,'DisplayName','Pref');
% plot(t,Q_res,'LineWidth',1.5,'DisplayName','Q');
% plot(t,Q1ref*ones(1,length(t)),'LineWidth',1.5,'DisplayName','Qref');
% % plot(t,abc_power,'LineWidth',1.5,'DisplayName','ABC')
% legend('Interpreter','latex')
% title('PQ Tracking')
% % ylim([-1000 11000])
% xlabel('Time (s)')
% ylabel('Power (W)')

% f4 = figure('Position',[200 100 550 400]);
% plot(t,vg_al,'DisplayName','$Vg_{alpha}$')
% hold on
% plot(t,vg_be,'DisplayName','$Vg_{beta}$')
% plot(t,vc_al,'DisplayName','$Vc_{alpha}$')
% plot(t,vc_be,'DisplayName','$Vc_{beta}$')
% xlim([0 0.05])
% ylim([-1000 1000])
% title('Alpha Beta Voltages')
% legend('Interpreter','latex')

f5 = figure('Position',[200 100 550 400]);
plot(t,P_res,'LineWidth',1.5,'DisplayName','P');
hold on
plot(t,Q_res,'LineWidth',1.5,'DisplayName','Q');
plot(t,Ptot_ref,'LineWidth',1.5,'DisplayName','$P_{ref}$')
plot(t,Qtot_ref,'LineWidth',1.5,'DisplayName','$Q_{ref}$')
% plot(t,abc_power,'LineWidth',1.5,'DisplayName','ABC')
legend('Interpreter','latex')
title('PQ Tracking')
% ylim([-1000 11000])
xlabel('Time (s)')
ylabel('Power (W)')


f6 = figure('Position',[200 100 550 400]);
plot(t,i_a)
hold on
plot(t,i_b)
xlabel('Time (s)')
ylabel('Current (A)')
title('Currents')
legend(["$I_a$","$I_b$"],'Interpreter','latex')

f7 = figure('Position',[200 100 550 400]);
% plot(t,vc_abc)
% hold on
plot(t,vg_abc)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Grid Phase Voltages')
legend(["$Vg_a$","$Vg_b$","$Vg_c$"],'Interpreter','latex')

f8 = figure('Position',[200 100 550 400]);
plot(t,res_eng)
xlabel('Time (s)')
ylabel('Energy (J)')
title('Transformed Capacitor Energy')
legend(["$E_{tot}$","$E_{a-b}$","$E_{a-c}$"],'Interpreter','latex')

f9 = figure('Position',[200 100 550 400]);
plot(t,e_eng)
xlabel('Time (s)')
ylabel('Energy (J)')
title('Energy Errors')
legend(["$e_{tot}$","$e_{a-b}$","$e_{a-c}$"],'Interpreter','latex')

% f10 = figure('Position',[200 100 550 400]);
% yyaxis right
% plot(t,alph-w*t,'LineWidth',1.5,'DisplayName','$\phi_{v_\alpha}$');
% hold on
% plot(t,beph-w*t,'LineWidth',1.5,'DisplayName','$\phi_{v_\beta}$');
% ylabel('Phase (rad)')
% yyaxis left
% plot(t,almag,'LineWidth',1.5,'DisplayName','$|v_\alpha|$');
% hold on
% plot(t,bemag,'LineWidth',1.5,'DisplayName','$|v_\beta|$');
% ylabel('Voltage (V)')
% xlabel('Time (s)')
% title('SOGI Tracking')
% legend(Interpreter='latex')

