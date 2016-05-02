clear all
close all



A1 = load('waveElevationY082.csv');
B1 = load('waveCutsOlivieri0082.csv');
figure(1)
hold off
plot(A1(:,20)/5.72-0.04+0.5,A1(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B1(:,1),B1(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.082','Fontsize',18)
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)
print -dpng -color waveElevationY082.png

A2 = load('waveElevationY172.csv');
B2 = load('waveCutsOlivieri0172.csv');
figure(2)
hold off
plot(A2(:,20)/5.72-0.04+0.5,A2(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B2(:,1),B2(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.172','Fontsize',18)
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)
print -dpng -color waveElevationY172.png

A3 = load('waveElevationY259.csv');
B3 = load('waveCutsOlivieri0259.csv');
figure(3)
hold off
plot(A3(:,20)/5.72-0.04+0.5,A3(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B3(:,1),B3(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.259','Fontsize',18)
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)
print -dpng -color waveElevationY259.png



A4 = load('waveElevationY347.csv');
B4 = load('waveCutsOlivieri0347.csv');
figure(4)
hold off
plot(A4(:,20)/5.72-0.04+0.5,A4(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B4(:,1),B4(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.347','Fontsize',18)
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)
print -dpng -color waveElevationY347.png

figure(5)
C = load('DTMB_Fr028_fine_test_water_line.txt');D1 = load('goteborgExpFr028_2340.csv');
D1 = load('goteborgExpFr028_2340.csv');
D2 = load('goteborgExpFr028_5415.csv');
D3 = load('goteborgExpFr028_5512.csv');
plot(C(:,1)/5.72+0.5-0.035,C(:,3)/5.72,'*')
hold on; grid on; plot(D1(:,1),D1(:,2)-((0.0158-(0.0158-0.005)*D1(:,1))/5.72),'r','Linewidth',2)
hold on; grid on; plot(D2(:,1),D2(:,2)-((0.0158-(0.0158-0.005)*D2(:,1))/5.72),'g','Linewidth',2)
hold on; grid on; plot(D3(:,1),D3(:,2)-((0.0158-(0.0158-0.005)*D3(:,1))/5.72),'m','Linewidth',2)
legend('WaveBEM','Exp. 2340', 'Exp. 5415', 'Exp. 5512')
title('DTMB-5415 Water Line','Fontsize',18)
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)
print -dpng -color waveElevationHull.png


figure(6)
subplot(4,1,1)
hold off
plot(A1(:,20)/5.72-0.04+0.5,A1(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B1(:,1),B1(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.082','Fontsize',14, 'Fontweight','Bold')
%xlabel('x/L_{pp}','Fontsize',10)
ylabel('z/L_{pp}','Fontsize',14)

subplot(4,1,2)
hold off
plot(A2(:,20)/5.72-0.04+0.5,A2(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B2(:,1),B2(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.172','Fontsize',14, 'Fontweight','Bold')
%xlabel('x/L_{pp}','Fontsize',10)
ylabel('z/L_{pp}','Fontsize',14)

subplot(4,1,3)
hold off
plot(A3(:,20)/5.72-0.04+0.5,A3(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B3(:,1),B3(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.259','Fontsize',14, 'Fontweight','Bold')
%xlabel('x/L_{pp}','Fontsize',10)
ylabel('z/L_{pp}','Fontsize',14)

subplot(4,1,4)
hold off
plot(A4(:,20)/5.72-0.04+0.5,A4(:,22)/5.72,'*')
axis([(-0.1*5.72)/5.72 (2.5*5.72)/5.72 -0.008 0.008])
hold on; grid on; plot(B4(:,1),B4(:,2),'r','Linewidth',2)
legend('WaveBEM','Experimental')
title('y/L_{pp} = 0.347','Fontsize',14, 'Fontweight','Bold')
xlabel('x/L_{pp}','Fontsize',14)
ylabel('z/L_{pp}','Fontsize',14)

print -dpng -color waveElevationDTMB4Cuts.png
