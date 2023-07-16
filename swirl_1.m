% Swirl analysis using tangential holes
clc;clear all;
close all;

set(groot,'defaultLineLineWidth',2.0)

vary = input("What do you want to vary? Type 1 for offset, 2 for orf rad, 3 for exit rad, 4 for num of orf, 5 for cd orifice: ");
start=input("Enter the starting value of parameter you want to vary in metres: ");
final=input("Enter the ending value of parameter you want to vary in metres: ");
data_points=input("Enter number of data points you want between them: ");

% thrust chamber parameter
chamberP=35; % bar

% target
target_deltaP=0.25*chamberP; % bar

% Injector element geometric parameters
Rbx=linspace(3.4e-3,3.4e-3,data_points); % offset
rbx=linspace(1.5e-3,1.5e-3,data_points); % orf rad
rc=linspace(5.5e-3,5.5e-3,data_points); % exit rad
n=linspace(6,6,data_points); % num of orf
cd_orf=linspace(0.65,0.65,data_points); % orifice cd

% flow rate
mdot_feed=2; %kg/s
mdot_per_element=mdot_feed./n;

% oxidizer property (considering liquid injection)
temp=273;
rho=linspace(917.57,917.57,data_points);
mu=linspace(90.442e-6,90.442e-6,data_points);

if vary==1
    Rbx=linspace(start,final,data_points);
    variable=Rbx;
elseif vary==2
    rbx=linspace(start, final, data_points);
    variable=rbx;
elseif vary==3
    rc=linspace(start, final, data_points);
    variable=rc;
elseif vary==4
    n=linspace(start, final, data_points);
    variable=n; 
elseif vary==5
    cd_orf=linspace(start, final, data_points);
    variable=cd_orf;
end

A = Rbx.*rc./(rbx.*rbx.*n);
orf_area=3.14.*rbx.*rbx;

for j=1:length(A)
fun = @(fio)func(fio,A(j)); % calling fsolve function to solve for fi
fi0=0.7; % initial fi guess
fio=fsolve(fun,fi0);
fo(j)=fio;
end

cd=sqrt(fo.^3./(2-fo));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fo)).^2-4.*cd.*cd.*A.*A))); % full spray cone angle
delp1=10^(-5).*(mdot_per_element./(cd_orf.*n.*orf_area.*sqrt(2*rho))).^2; % due to orifice
delp2=10^(-5).*(mdot_per_element./(cd.*3.14.*rc.*rc.*sqrt(2.*rho))).^2; % due to swirl with effects of friction included in cd
delp_tot=delp2+delp1;
axial_vel=mdot_per_element./(rho.*3.14.*rc.*rc.*fo);
tang_vel=axial_vel.*tand(alpha./2);
resultant_vel=sqrt(tang_vel.^2+tang_vel.^2);
calc=[cd,alpha, delp1, delp2, delp_tot,axial_vel,tang_vel,resultant_vel];

figure(1)
plot(A,cd)
hold on
scatter(A,cd)
hold on
plot(A,alpha.*0.01)
hold on
scatter(A,alpha.*0.01)
hold on
plot(A,fo)
hold on
scatter(A,fo)
grid on
xlabel("Swirl Number")
ylabel("Cd,fi,alpha")

figure(2)
plot(A,delp1)
hold on
scatter(A,delp1)
grid on
xlabel("Swirl Number")
ylabel("deltaP1-orf (bar)")

figure(3)
plot(A,delp2)
hold on
scatter(A,delp2)
grid on
xlabel("Swirl Number")
ylabel("deltaP2-due to swirl (bar)")

figure(4)
plot(variable, A)
hold on
scatter(variable, A)
grid on
ylabel("Swirl Number")
xlabel("Variable")

figure(5)
plot(A,delp_tot)
hold on
scatter(A,delp_tot)
grid on
xlabel("Swirl Number")
ylabel("deltaP_total (bar)")

figure(6)
plot(A,axial_vel)
hold on
scatter(A,axial_vel)
grid on
xlabel("Swirl Number")
ylabel("axial velocity (m/s)")

figure(7)
plot(A,tang_vel)
hold on
scatter(A,tang_vel)
grid on
xlabel("Swirl Number")
ylabel("tangential velocity (m/s)")

figure(8)
plot(A,resultant_vel)
hold on
scatter(A,resultant_vel)
grid on
xlabel("Swirl Number")
ylabel("resultant velocity (m/s)")


% function
function F = func(f,A)
F=1-sqrt((A.^2./(1-f))+(f.^(-2))).*sqrt(f.^3/(2-f));
end



