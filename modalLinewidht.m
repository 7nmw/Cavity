clc
clear all
close all
%1 Khz for 1 mode
% 
ni=2.2111;%refractive index idler photon
ns=2.2540; %signal photon
L=20*10^(-3); % m crystall size
c0=299792458; % m/sec
% Dphy=49*10^(-3); % m distans cavity
Dphy=[11:1:251]*10^(-3); % m distans cavity
tp=3*10^(-2); % коэффициент пропуская первого зеркала (в отн.единицах)
tout=3*10^(-2); % коэффициент пропуская 2-ого зеркала
alfa=7*10^(-2); % Измеренные внутренние потери резонатора на длине волны резонирующего сигнала
nair=1;
r=120*10^(-3) %  радиус кривизны зеркал в метрах
lambda = 867*10^(-9); % метра

%
% 
d=Dphy-L+L/ns; % metr effective distance
dcry=Dphy+L-L/ns;

r0=(L/c0)*(ni-ns);% время жизни фотона в резонатора
Twpeaks=abs(r0*10^12) %temporal linewidth ps

Tr=(nair*2*(dcry-L)+ns*2*L)/c0; % время полного обхода резонатора туда и обратно в секундах
Tr1=Tr*10^9; % в наносекундах

alfa_rt=tp+tout+2*alfa; %The power loss from the cavity during a round-trip
alfaCav=alfa_rt./Tr; %The rate at which photons are lost from the cavit
% delta_Omega_mode==delta_Omega_cav==alfaCav; %angular frequency width of the signal and idler modes
Delta_V_mode=alfaCav./2*pi; %linewidth


Delta_V_FSR=1./Tr; % межмодовое расстояние в Hz/ The cavity free spectral range


%Finesse1=2*pi/(tp+tout+2*alfa) % Finesse
Finesse=2*pi/(1-(0.97*0.97*(1-2*0.07))) % Finesse


Vmode=Delta_V_FSR/Finesse; %ШИрина модыMHz

%%%%%%%%%%%%%%%%


% d1=d*10^3

U=1-2*(d/r)+(d/r).^2; %stability region of the cavity.

wcf=sqrt((lambda*r)/(2*pi)); %focused beam waist radius in a confocal symmetric
w0=sqrt(sqrt(((d/r).*(1-d/r).*(2-d/r)*wcf)./(1-d/r))); %The minimum waist, 𝑤𝑤0, located in the center of the cavity is given by
d1=d*10^3;
w0=abs(w0);
w0=w0*10^3;
wop=sqrt(L*lambda/2*pi*2.84)*10^6; %pump waist radius (according to [54]) for optimal interaction inside a bulk crystal
% plot(d1,U)
dcry = dcry*10^3;
Vmode=Vmode.*10^(-6);
Delta_V_FSR=Delta_V_FSR.*10^(-9);


% Top plot
ax1 = nexttile;
plot(dcry, U)
ylabel('U')
xlabel('d_{cry}, мм')
grid on
% legend([h1,h2,h3,h4],{'Data','Best w/ rMax = 30'},'Location','SouthEast');
% text(15,1,sprintf('Finesse = %5.2f \n ReflectionOut = %5.2f \n ReflectionIn = %5.3f \n Losses = %5.2f \n temporal linewidth ps = %5.2f',Finesse,tout,tp,alfa,Twpeaks),'HorizontalAlignment','left','FontSize',12);
% 
% Bottom plot
ax2 = nexttile;
plot(dcry,w0)
title(ax2,'focused beam waist')
ylabel(ax2,'w0 (um)')
xlabel(ax2,'d_{cry} (mm)')
grid on
% axis( [ 0, 51, 0, 100 ] )
% 
% 
% Bottom plot
ax3 = nexttile;
plot(dcry,Vmode)
title(ax3,'modal linewidth')
ylabel(ax3,'Vmode (MHz)')
xlabel(ax3,'d_{cry} (mm)')
grid on
% text(10,35,sprintf('Finesse = %5.2f \n ReflectionOut = %5.2f \n ReflectionIn = %5.3f \n Losses = %5.2f \n temporal linewidth ps = %5.2f',Finesse,tout,tp,alfa,Twpeaks),'HorizontalAlignment','left','FontSize',12);
% 


% Bottom plot
ax4 = nexttile;
plot(dcry,Delta_V_FSR)
title(ax4,'DeltaFSR')
ylabel(ax4,'FSR (GHz)')
xlabel(ax4,'d_{cry} (mm)')
grid on


% SPF=(sinc(pi*r0.*Vmode)).^2;
% G=exp(-alfa*
% delay=r/Tr
% figure
% plot(SPF)
%delta_w=0.88*pi/r0;

%delta_V=delta_w/(2*pi);



