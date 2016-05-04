close all
clear all

L = 2.5;
Sw = 0.933845;
rho = 1025.1;
g = 9.81;
Fr = [0.250 0.267 0.289 0.316 0.354 0.408];
v = Fr*sqrt(g*L);


exp250 = load('./tokyoFr025.csv');
num250 = load('./test_wigley_Fr250_water_line.txt');
[S,I] = sort(num250(:,1));

num250Sorted = num250(I,:);
figure(1)
plot(exp250(:,1),exp250(:,2),'k*')
hold on
plot(num250Sorted(:,1)/L,2*g*num250Sorted(:,3)/v(1)^2,'r')
grid on

exp267 = load('./tokyoFr0267.csv');
num267 = load('./test_wigley_Fr267_water_line.txt');
[S,I] = sort(num267(:,1));

num267Sorted = num267(I,:);
figure(2)
plot(exp267(:,1),exp267(:,2),'k*')
hold on
plot(num267Sorted(:,1)/L,2*g*num267Sorted(:,3)/v(2)^2,'r')
grid on


exp289 = load('./tokyoFr0289.csv');
num289 = load('./test_wigley_Fr289_water_line.txt');
[S,I] = sort(num289(:,1));

num289Sorted = num289(I,:);
figure(3)
plot(exp289(:,1),exp289(:,2),'k*')
hold on
plot(num289Sorted(:,1)/L,2*g*num289Sorted(:,3)/v(3)^2,'r')
grid on

exp316 = load('./tokyoFr0316.csv');
num316 = load('./test_wigley_Fr316_water_line.txt');
[S,I] = sort(num316(:,1));

num316Sorted = num316(I,:);
figure(4)
plot(exp316(:,1),exp316(:,2),'k*')
hold on
plot(num316Sorted(:,1)/L,2*g*num316Sorted(:,3)/v(4)^2,'r')
grid on

exp354 = load('./tokyoFr0354.csv');
num354 = load('./test_wigley_Fr354_water_line.txt');
[S,I] = sort(num354(:,1));

num354Sorted = num354(I,:);
figure(5)
plot(exp354(:,1),exp354(:,2),'k*')
hold on
plot(num354Sorted(:,1)/L,2*g*num354Sorted(:,3)/v(5)^2,'r')
grid on



exp408 = load('./tokyoFr0408.csv');
num408 = load('./test_wigley_Fr408_water_line.txt');
[S,I] = sort(num408(:,1));

num408Sorted = num408(I,:);
figure(6)
plot(exp408(:,1),exp408(:,2),'k*')
hold on
plot(num408Sorted(:,1)/L,2*g*num408Sorted(:,3)/v(6)^2,'r')
grid on


forces250 = load('./test_wigley_Fr250_force.txt');
forces267 = load('./test_wigley_Fr267_force.txt');
forces289 = load('./test_wigley_Fr289_force.txt');
forces316 = load('./test_wigley_Fr316_force.txt');
forces354 = load('./test_wigley_Fr354_force.txt');
forces408 = load('./test_wigley_Fr408_force.txt');

figure(7)
plot(forces250(:,1),forces250(:,2),forces267(:,1),forces267(:,2),forces289(:,1),forces289(:,2),forces316(:,1),forces316(:,2),forces354(:,1),forces354(:,2),forces408(:,1),forces408(:,2))
grid on

Dw = [mean(forces250(end-30:end,2)),
      mean(forces267(end-30:end,2)),
      mean(forces289(end-30:end,2)),
      mean(forces316(end-30:end,2)),
      mean(forces354(end-30:end,2)),
      mean(forces408(end-30:end,2))]

Cdw = Dw./(0.5*rho*v.^2*Sw)'

CdwExp1 = load('wigleyDragExp1.csv');
CdwExp2 = load('wigleyDragExp2.csv');
CdwExp3 = load('wigleyDragExp3.csv');
CdwExp4 = load('wigleyDragExp4.csv');

figure(8)
plot(Fr,Cdw' ,'r*-')
hold on

plot(CdwExp1(:,1),CdwExp1(:,2),'k*')
plot(CdwExp2(:,1),CdwExp2(:,2),'kd')
plot(CdwExp3(:,1),CdwExp3(:,2),'ko')
plot(CdwExp4(:,1),CdwExp4(:,2),'k+')
axis([0.08 0.42 -1e-3 3.8e-3])
grid on

