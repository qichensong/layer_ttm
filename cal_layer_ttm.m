close all
T = 295;
C = getC('Al',T);
kz = 10;



wmin = 5e3;
wmax = 5e9;
nw = 40;
omega = exp(linspace(log(wmin),log(wmax),nw))*2*pi;

w0 = 20e-6;
w1 = 20e-6;

alp = 1/(16e-9);
gamma = 1/(17e-9);

Ce = 3e4;
kez = getk('Al',T)-10;
g = 5e16;

L = [100e-9 1e-3];
Gep = 1e2;
Gpp = 3e8;
k2x = getk('Si',T);
k2z = getk('Si',T);
C2 = getC('Si',T);
%g = 5e12;


data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
amp = reshape(data(1,1,:,1),[length(omega) 1]);
figure(1)
hold on
plot(omega/2/pi,amp,'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
figure(2)
hold on
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')

tdelay = linspace(500e-12,8e-9,100);
datat = return_tdtr_real_space_layer_ttm(tdelay,6e6*2*pi,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
figure(3)
amp = reshape(datat(1,1,:,1),[length(tdelay) 1]);
hold on
plot(tdelay,amp,'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Delay time (s)')
ylabel('Amplitude')
figure(4)
hold on
plot(tdelay,reshape(datat(1,1,:,2),[length(tdelay) 1]),'LineStyle','--','DisplayName','w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Delay time (s)')
ylabel('Phase (deg)')





Gep1 = 3e8;


data = return_data_real_space_layer_ttm(omega,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
amp = reshape(data(1,1,:,1),[length(omega) 1]);
figure(1)
hold on
plot(omega/2/pi,amp,'LineStyle','-','DisplayName','with G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Frequency (Hz)')
figure(2)
hold on
plot(omega/2/pi,reshape(data(1,1,:,2),[length(omega) 1]),'LineStyle','-','DisplayName','with G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Frequency (Hz)')

tdelay = linspace(500e-12,8e-9,100);
datat = return_tdtr_real_space_layer_ttm(tdelay,6e6*2*pi,w0,w1,[0],[0],C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
figure(3)
amp = reshape(datat(1,1,:,1),[length(tdelay) 1]);
hold on
plot(tdelay,amp,'LineStyle','-','DisplayName','with w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Delay time (s)')
figure(4)
hold on
plot(tdelay,reshape(datat(1,1,:,2),[length(tdelay) 1]),'LineStyle','-','DisplayName','with w/o G_{ep}','LineWidth',1.5)
set(gca, 'XScale', 'log')
xlabel('Delay time (s)')
for i = 1:4
    figure(i)
    box on
    legend boxoff
    %saveas(gcf,[num2str(i) '.svg'])
end


zz = linspace(0,100e-9,300);
zz2 = linspace(0,100e-9,200);
[te,tp,t2]=return_data_real_space_z_layer_ttm(omega,w0,w1,[0],[0],zz,zz2,C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep,Gpp,k2x,k2z,C2);
figure(5)
iw = 20;
% plot(zz,te(:,iw,1).*sin(te(:,iw,2)/180*pi))
hold on
% plot(zz,tp(:,iw,1).*sin(tp(:,iw,2)/180*pi))
% plot(zz(end)+zz2,t2(:,iw,1).*sin(t2(:,iw,2)/180*pi))
plot(zz,te(:,iw,1).*cos(te(:,iw,2)/180*pi),'r--','LineWidth',2)
plot(zz,tp(:,iw,1).*cos(tp(:,iw,2)/180*pi),'g--','LineWidth',2)
plot(zz(end)+zz2,t2(:,iw,1).*cos(t2(:,iw,2)/180*pi),'b--','LineWidth',2)
figure(6)
plot(zz,te(:,iw,1),'r--','LineWidth',2)
hold on
plot(zz,tp(:,iw,1),'g--','LineWidth',2)
plot(zz(end)+zz2,t2(:,iw,1),'b--','LineWidth',2)
[te,tp,t2]=return_data_real_space_z_layer_ttm(omega,w0,w1,[0],[0],zz,zz2,C,kz,kz,alp,gamma,Ce,kez,kez,g,L,Gep1,Gpp,k2x,k2z,C2);
figure(5)
% plot(zz,te(:,iw,1).*sin(te(:,iw,2)/180*pi))
hold on
% plot(zz,tp(:,iw,1).*sin(tp(:,iw,2)/180*pi))
% plot(zz(end)+zz2,t2(:,iw,1).*sin(t2(:,iw,2)/180*pi))
plot(zz,te(:,iw,1).*cos(te(:,iw,2)/180*pi),'r-','LineWidth',2)
plot(zz,tp(:,iw,1).*cos(tp(:,iw,2)/180*pi),'g-','LineWidth',2)
plot(zz(end)+zz2,t2(:,iw,1).*cos(t2(:,iw,2)/180*pi),'b-','LineWidth',2)
figure(6)
plot(zz,te(:,iw,1),'r-','LineWidth',2)
hold on
plot(zz,tp(:,iw,1),'g-','LineWidth',2)
plot(zz(end)+zz2,t2(:,iw,1),'b-','LineWidth',2)