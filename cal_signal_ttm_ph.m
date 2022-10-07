close all
T = 300;
C = getC('Ge',T);
kz = getk('Ge',T);

wmin = 5e3;
wmax = 5e8;
omega = exp(linspace(log(wmin),log(wmax),100))*2*pi;

w0 = 10e-6;
w1 = 10e-6;

alp = 1/(16e-9);
gamma = 1/(17e-9);

Ce = 3e4;
kez = 0.1;
g = 5e16;
%g = 5e12;


data = return_data_real_space_ttm(omega,w0,w1,[0],[0],C,kz,0.1*kz,alp,gamma,Ce,kez,kez,g);
figure(1)
hold on
plot(omega/2/pi,reshape(data(1,1,:,1),[length(omega) 1])/data(1,1,1,1),'LineStyle','--','DisplayName','el')
set(gca, 'XScale', 'log')
figure(2)
hold on
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','--','DisplayName','el')
set(gca, 'XScale', 'log')

data = return_data_real_space_ttm_ph(omega,w0,w1,[0],[0],C,kz,0.1*kz,alp,gamma,Ce,kez,kez,g);
figure(1)
plot(omega/2/pi,reshape(data(1,1,:,1),[length(omega) 1])/data(1,1,1,1),'LineStyle','-','DisplayName','ph')
figure(2)
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','-','DisplayName','ph')

data = return_data_real_space_ttm_ph(omega,w0,w1,[0],[0],C,kz,0.1*kz,alp,gamma,Ce,kez,kez,g*1e8);
figure(1)
plot(omega/2/pi,reshape(data(1,1,:,1),[length(omega) 1])/data(1,1,1,1),'LineStyle','-.','DisplayName','Fourier')
legend
figure(2)
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','-.','DisplayName','Fourier')
legend
