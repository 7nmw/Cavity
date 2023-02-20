clc
close all
clear all
x=[-4:0.02:4]*1e-3; 
y=[-4:0.02:4]*1e-3;  
[X Y]=meshgrid(x,y);
[phi,ro] = cart2pol(X,Y);

w0=0.001; % beam spot m
lyam=532e-9; % ������ ����� m; 
k=(2*pi)/lyam; %wavenumber; 


%��� ������������� ������, ������� 1, ����� ���� �� ����� P,Q ������
% P=-2;
% Q=5;
% m=9; % ����������� ������
% ka1=m/(2*Q);
% ka1=fix(ka1);
% N1=(P+Q)*ka1;
% N2=m-2*Q*ka1
% 
% ka2=m/abs(P-Q);
% ka2=fix(ka2);
% N3=m+(P-Q)*ka2;
% N4=m-2*Q*ka2

%��� ������������� ������, �������� 2, ����� ��� ����� P, Q �������� 
% P=1;
% Q=5;
% m=2; % ����������� ������
% ka1=m/(Q)
% ka1=fix(ka1)
% N1=(P+Q)*ka1;
% N2=m-Q*ka1
% 
% ka2=m/abs(P-Q);
% ka2=fix(ka2);
% N3=m+(P-Q)*ka2;
% N4=m-Q*ka2

N2=-1
N4=5

%���������� �������, ���������� � ������������ �������
% L=2;
p=0;
p2=0;
R=0.003; %metr
l=0.000508; %���������� ��������� M1 and M2, M1 and M3 � ������
lo=0.001; %����� ����� ��������� 
Lo=l*2+lo; %������ ������� ������ ����������
tetao=1/5; %�������� �������� �����
J=2;

A=1-4*l/R-2*lo/R+4*l*lo/R^2;
D=1-4*l/R-2*lo/R+4*l*lo/R^2;
teta=tetao*acos(A+D/2); % ��� ������ ���������� ����� �������������� �� ����
gamma=2*p+abs(J)+2*teta*J+1;

% syms m
% PolinomL=symsum((-1)^m.*((factorial(L1+p1))./((factorial(p1-m))*(factorial(L1+m))*(factorial(m))))*(2*ro.^2/w0^2).^m, m, 0, p1); %������� ������� ������� (2)
P=(exp(-ro.^2/(w0^2))).*(exp(-1i*N2.*phi)); 
Lg=(sqrt(2*factorial(p)/(pi*(factorial(abs(N2)+p)))))*(1/w0)*((sqrt(2)*ro./w0).^abs(N2)).*P; %���� ������-������
Lg=double(Lg);


syms m
PolinomL=symsum((-1)^m.*((factorial(abs(N4)+p2))./((factorial(p2-m))*(factorial(abs(N4)+m))*(factorial(m))))*(2*ro.^2/w0^2).^m, m, 0, p2); %������� ������� ������� (2)
P1=(exp(-ro.^2/(w0^2))).*(exp(-1i*N4.*phi)); 
Lg1=(sqrt(2*factorial(p2)/(pi*(factorial(abs(N4)+p2)))))*(1/w0)*((sqrt(2)*ro./w0).^abs(N4)).*P1.*PolinomL; %���� ������-������
Lg1=double(Lg1);

% figure
% mesh(X,Y,abs(Lg));
% figure
% mesh(X,Y,abs(Lg1));

Lg2=Lg+Lg1;
figure
mesh(X,Y,abs(Lg2));


% phaseshift=exp(1i*k*Lo-1i*gamma*acos(A+D/2)); % ����� ����
% beamangle=exp(teta*acos(A+D/2)); %

% syms L
% Fo=symsum(((sqrt(2*factorial(p)/(pi*(factorial(L+p)))))*(1/w0)*((sqrt(2)*ro./w0).^L)).*((exp(-ro.^2/(w0^2))).*(exp(-1i*L.*phi))), L, 0, 2);

% Fo=symsum((exp(-ro.^2).*(ro.^L).*(exp(-1i*L.*phi))), L, 0, 1);

% Fo=double(Fo);
% Fo=char(Fo)
% F=phaseshift.*beamangle.*Fo.*exp(1i*J*teta);

% s=class(phaseshift)
% F=double(F); %Convert Symbolic Number to Double Precision
% C3=char(C)
% mesh(X,Y,abs(F));
% set(gca,'view',[0 90]);grid on;