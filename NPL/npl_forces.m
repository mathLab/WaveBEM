clear all
close all
Fr = [0.2 0.25 0.3 0.35 0.4];

Sw = 149.611;
Lw = 30.5002;
g = 9.81;
rho = 1025;

V = Fr*sqrt(Lw*g)

experiments = load('npl_exp_drag.csv');

CR_exp = experiments(:,2);
Fr_exp = experiments(:,1);
Fr_exp = Fr_exp*0.5144444/sqrt(0.3048*9.81);

CR = interp1(Fr_exp,CR_exp,Fr,'cubic','extrap')


figure(2)
plot(Fr_exp,CR_exp,"k--","linewidth",4)
hold on
%plot(Fr,CR,'b*-')
grid on

figure(1)
R_exp = 1/2*Sw*rho*V.^2.*CR
plot(Fr,R_exp,"ko--","linewidth",4)
grid on

%%%%%%%%%%%%%%%% Fr = 0.2
WaveBEM_Fr02_old = load('npl_Fr02_force.txt');
WaveBEM_Fr02_new = load('npl_Fr020_damped_less_force.txt');

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr02_old)
    if (WaveBEM_Fr02_old(ii,1) > 40)
       tmp(count) = WaveBEM_Fr02_old(ii,2)-WaveBEM_Fr02_old(ii,5)+WaveBEM_Fr02_old(ii,8)+WaveBEM_Fr02_old(ii,11);
    endif
end
drag_WaveBEM_Fr02_old = mean(tmp)

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr02_new)
    if (WaveBEM_Fr02_new(ii,1) > 60 && WaveBEM_Fr02_new(ii,1) < 90)
       tmp(count) = WaveBEM_Fr02_new(ii,2)-WaveBEM_Fr02_new(ii,5)+WaveBEM_Fr02_new(ii,8)+WaveBEM_Fr02_new(ii,11);
    endif
end
drag_WaveBEM_Fr02_new = mean(tmp)

%%%%%%%%%%%%%%%% 



%%%%%%%%%%%%%%%% Fr = 0.25
WaveBEM_Fr025_old = load('npl_new_mesh_test_force.txt');
WaveBEM_Fr025_new = load('npl_Fr025_damped_less_new_force.txt');

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr025_old)
    if (WaveBEM_Fr025_old(ii,1) > 40)
       tmp(count) = WaveBEM_Fr025_old(ii,2)-WaveBEM_Fr025_old(ii,5)+WaveBEM_Fr025_old(ii,8)+WaveBEM_Fr025_old(ii,11);
    endif
end
drag_WaveBEM_Fr025_old = mean(tmp)

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr025_new)
    if (WaveBEM_Fr025_new(ii,1) > 60 && WaveBEM_Fr025_new(ii,1) < 90)
       tmp(count) = WaveBEM_Fr025_new(ii,2)-WaveBEM_Fr025_new(ii,5)+WaveBEM_Fr025_new(ii,8)+WaveBEM_Fr025_new(ii,11);
    endif
end
drag_WaveBEM_Fr025_new = mean(tmp)

%%%%%%%%%%%%%%%% 



%%%%%%%%%%%%%%%% Fr = 0.3
WaveBEM_Fr03_old = load('npl_Fr03_new_force.txt');
WaveBEM_Fr03_new = load('npl_Fr030_damped_less_force.txt');

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr03_old)
    if (WaveBEM_Fr03_old(ii,1) > 60 && WaveBEM_Fr03_old(ii,1) < 90)
       tmp(count) = WaveBEM_Fr03_old(ii,2)-WaveBEM_Fr03_old(ii,5)+WaveBEM_Fr03_old(ii,8)+WaveBEM_Fr03_old(ii,11);
    endif
end
drag_WaveBEM_Fr03_old = mean(tmp)

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr03_new)
    if (WaveBEM_Fr03_new(ii,1) > 60 && WaveBEM_Fr03_new(ii,1) < 90)
       tmp(count) = WaveBEM_Fr03_new(ii,2)-WaveBEM_Fr03_new(ii,5)+WaveBEM_Fr03_new(ii,8)+WaveBEM_Fr03_new(ii,11);
    endif
end
drag_WaveBEM_Fr03_new = mean(tmp)

%%%%%%%%%%%%%%%% 




%%%%%%%%%%%%%%%% Fr = 0.35
WaveBEM_Fr035_old = load('npl_Fr035_new_force.txt');
WaveBEM_Fr035_new = load('npl_Fr035_damped_less_force.txt');

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr035_old)
    if (WaveBEM_Fr035_old(ii,1) > 40)
       tmp(count) = WaveBEM_Fr035_old(ii,2)-WaveBEM_Fr035_old(ii,5)+WaveBEM_Fr035_old(ii,8)+WaveBEM_Fr035_old(ii,11);
    endif
end
drag_WaveBEM_Fr035_old = mean(tmp)

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr035_new)
    if (WaveBEM_Fr035_new(ii,1) > 60 && WaveBEM_Fr035_new(ii,1) < 90)
       tmp(count) = WaveBEM_Fr035_new(ii,2)-WaveBEM_Fr035_new(ii,5)+WaveBEM_Fr035_new(ii,8)+WaveBEM_Fr035_new(ii,11);
    endif
end
drag_WaveBEM_Fr035_new = mean(tmp)

%%%%%%%%%%%%%%%% 




%%%%%%%%%%%%%%%% Fr = 0.4
WaveBEM_Fr04_old = load('npl_Fr04_damped_force.txt');
WaveBEM_Fr04_new = load('npl_Fr04_damped_less_force.txt');

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr04_old)
    if (WaveBEM_Fr04_old(ii,1) > 60 && WaveBEM_Fr04_old(ii,1) < 90)
       tmp(count) = WaveBEM_Fr04_old(ii,2)-WaveBEM_Fr04_old(ii,5)+WaveBEM_Fr04_old(ii,8)+WaveBEM_Fr04_old(ii,11);
    endif
end
drag_WaveBEM_Fr04_old = mean(tmp)

tmp = [];
count = 1;
for ii=1:length(WaveBEM_Fr04_new)
    if (WaveBEM_Fr04_new(ii,1) > 60 && WaveBEM_Fr04_new(ii,1) < 90)
       tmp(count) = WaveBEM_Fr04_new(ii,2)-WaveBEM_Fr04_new(ii,5)+WaveBEM_Fr04_new(ii,8)+WaveBEM_Fr04_new(ii,11);
    endif
end
drag_WaveBEM_Fr04_new = mean(tmp)

%%%%%%%%%%%%%%%% 

numerical_Doctors = [drag_WaveBEM_Fr02_old drag_WaveBEM_Fr025_old drag_WaveBEM_Fr03_old drag_WaveBEM_Fr035_old];
numerical_Doctors_mod = [drag_WaveBEM_Fr02_new drag_WaveBEM_Fr025_new drag_WaveBEM_Fr03_new drag_WaveBEM_Fr035_new drag_WaveBEM_Fr04_new];


figure(1)
hold on
%plot(Fr(1:4),numerical_Doctors,'k*-')
plot(Fr,numerical_Doctors_mod,"kd-","linewidth",4, "MarkerFaceColor", "k")
title('NPL Hull Total Drag','Fontsize',24)
xlabel('Fr','Fontsize',18)
ylabel('D [N]','Fontsize',18)
legend('Experimental','WaveBEM')
print('NPL_drag_curve.png','-dpng')

figure(2)
%plot(Fr(1:4),numerical_Doctors./(1/2*Sw*rho*V(1:4).^2),'k*-')
plot(Fr,numerical_Doctors_mod./(1/2*Sw*rho*V.^2),"kd-","linewidth",4, "MarkerFaceColor", "k")
title('NPL Hull Total Drag Coefficient','Fontsize',24)
xlabel('Fr','Fontsize',18)
ylabel('C_D [N]','Fontsize',18)
legend('Experimental','WaveBEM')
print('NPL_drag_coeff_curve.png','-dpng')


force_errors = abs(R_exp-numerical_Doctors_mod)./R_exp


Cd_errors = abs(CR-numerical_Doctors_mod./(1/2*Sw*rho*V.^2))./CR
