% clear
% clc
% close all

k = sqrt(2);
w = 2*pi*50;
s = tf('s');
Q = 3;

alph = k*w*s/(s^2 + k*w*s + w^2);
beta = (k*w^2)/(s^2 + k*w*s + w^2);
% figure(1)
% bode(alph)
% figure(2)
% bode(beta)
%% N
s = tf('s');
N = (s^2 + w^2)/(s^2+2*s*w/Q+w^2);
figure(3)
bode(N)