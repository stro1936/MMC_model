clear
clc
close all 

syms ial ibe val vbe
% create clarke power expressions
ia = ial;
ib = (1/2)*(-ial + sqrt(3)*ibe);
ic = (1/2)*(-ial - sqrt(3)*ibe);

va = val;
vb = (1/2)*(-val + sqrt(3)*vbe);
vc = (1/2)*(-val - sqrt(3)*vbe);

Sa = expand(ia*va);
Sb = expand(ib*vb);
Sc = expand(ic*vc);   

albet_mat = [ial*val;ial*vbe;ibe*val;ibe*vbe]; % turn to matrix form

q = sqrt(3)/4;
M = [1 0 0 0; 1/4 -q -q 3/4; 1/4 q q 3/4;];

Sold = M*albet_mat;
C = [1 1 1;1 -1 0; 1 0 -1;]; 
Snew = C * Sold; %combination of powers we're interested in

syms ialmag ialph ibemag ibeph valmag valph vbemag vbeph ialx ialy ibex ibey valx valy vbex vbey real rational
syms Q1s P1s P2s P3s
% valmag = 565.69;
% valph = 1.572;
% vbemag = 565.69;
% vbeph = -3.143;
% Pr = 10e3;
% Qr = 0;
% P2r = 0;
% P3r = 0;

% P1, Q1, P2, & P3 are derived by hand from Snew from al be -> phase + magnitude
P1 = (3/2)*(valmag*ialmag/2)*cos(valph - ialph) + (3/2)*(vbemag*ibemag/2)*cos(vbeph - ibeph);
P1 = expand(P1);
P1 = subs(P1,ialmag*sin(ialph),ialy);
P1 = subs(P1,ibemag*sin(ibeph),ibey);
P1 = subs(P1,ialmag*cos(ialph),ialx);
P1 = subs(P1,ibemag*cos(ibeph),ibex);
P1 = subs(P1,valmag*sin(valph),valy);
P1 = subs(P1,vbemag*sin(vbeph),vbey);
P1 = subs(P1,valmag*cos(valph),valx);
P1 = subs(P1,vbemag*cos(vbeph),vbex);
P1 = expand(P1);
Q1 = (3/2)*(valmag*ialmag/2)*sin(valph - ialph) + (3/2)*(vbemag*ibemag/2)*sin(vbeph - ibeph);
Q1 = expand(Q1);
Q1 = subs(Q1,ialmag*sin(ialph),ialy);
Q1 = subs(Q1,ibemag*sin(ibeph),ibey);
Q1 = subs(Q1,ialmag*cos(ialph),ialx);
Q1 = subs(Q1,ibemag*cos(ibeph),ibex);
Q1 = subs(Q1,valmag*sin(valph),valy);
Q1 = subs(Q1,vbemag*sin(vbeph),vbey);
Q1 = subs(Q1,valmag*cos(valph),valx);
Q1 = subs(Q1,vbemag*cos(vbeph),vbex);
Q1 = expand(Q1);
P2 = (3/4)*(valmag*ialmag/2)*cos(valph - ialph) - (3/4)*(vbemag*ibemag/2)*cos(vbeph - ibeph) + ...
     (sqrt(3)/4)*(ialmag*vbemag/2)*cos(vbeph - ialph) + (sqrt(3)/4)*(ibemag*valmag/2)*cos(valph - ibeph);
P2 = expand(P2);
P2 = subs(P2,ialmag*sin(ialph),ialy);
P2 = subs(P2,ibemag*sin(ibeph),ibey);
P2 = subs(P2,ialmag*cos(ialph),ialx);
P2 = subs(P2,ibemag*cos(ibeph),ibex);
P2 = subs(P2,valmag*sin(valph),valy);
P2 = subs(P2,vbemag*sin(vbeph),vbey);
P2 = subs(P2,valmag*cos(valph),valx);
P2 = subs(P2,vbemag*cos(vbeph),vbex);
P2 = expand(P2);
P3 = (3/4)*(valmag*ialmag/2)*cos(valph - ialph) - (3/4)*(vbemag*ibemag/2)*cos(vbeph - ibeph) - ...
     (sqrt(3)/4)*(ialmag*vbemag/2)*cos(vbeph - ialph) - (sqrt(3)/4)*(ibemag*valmag/2)*cos(valph - ibeph);
P3 = expand(P3);
P3 = subs(P3,ialmag*sin(ialph),ialy);
P3 = subs(P3,ibemag*sin(ibeph),ibey);
P3 = subs(P3,ialmag*cos(ialph),ialx);
P3 = subs(P3,ibemag*cos(ibeph),ibex);
P3 = subs(P3,valmag*sin(valph),valy);
P3 = subs(P3,vbemag*sin(vbeph),vbey);
P3 = subs(P3,valmag*cos(valph),valx);
P3 = subs(P3,vbemag*cos(vbeph),vbex);
P3 = expand(P3);

eqs = [P1s==P1,Q1s==Q1,P2s==P2,P3s==P3]; % use P1 symbol to make matrix
vars = [ialx ialy ibex ibey];

[A,B] = equationsToMatrix(eqs,vars);
Ap1 = simplify(A\B) % these equations are implemented in Simulink model
Ap2 = simplify(inv(A)*B);
if simplify(Ap1-Ap2) ~= 0
    disp('inequality error')
end
%sol = solve(eqs,vars)
% ialmag_res= sol.ialmag
% ialph_res = sol.ialph
% ibemag_res = sol.ibemag
% ibeph_res = sol.ibeph
% F =@(x) [(3/2)*(valmag*x(1)/2)*cos(valph - x(2)) + (3/2)*(vbemag*x(3)/2)*cos(vbeph - x(4)) - Pr;
%          (3/2)*(valmag*x(1)/2)*sin(valph - x(2)) + (3/2)*(vbemag*x(3)/2)*sin(vbeph - x(4)) - Qr;
%          (3/4)*(valmag*x(1)/2)*cos(valph - x(2)) - (3/4)*(vbemag*x(3)/2)*cos(vbeph - x(4)) + ...
%          (sqrt(3)/4)*(x(1)*vbemag/2)*cos(vbeph - x(2)) + (sqrt(3)/4)*(x(3)*valmag/2)*cos(valph - x(4)) - P2r;
%      (3/4)*(valmag*x(1)/2)*cos(valph - x(2)) - (3/4)*(vbemag*x(3)/2)*cos(vbeph - x(4)) - ...
%      (sqrt(3)/4)*(x(1)*vbemag/2)*cos(vbeph - x(2)) - (sqrt(3)/4)*(x(3)*valmag/2)*cos(valph - x(4)) - P3r;];
% 
% x0 = [1;0;1;-pi/2];
% y = fsolve(F,x0,optimoptions('fsolve','Algorithm','levenberg-marquardt'))