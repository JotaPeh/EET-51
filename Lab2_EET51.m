clear; clc;
y=lfsr([1 0 0 0 0 0 1 1],[0 0 0 0 0 0 1]) %g(x)=xˆ5+xˆ2+1
N=127;%Period of the sequence N=2ˆL-1
Ryy=sequence_correlation(y,y,-130,130)%normalized autocorr.

p = plot(-130:1:130,Ryy);   %Plot autocorrelation
grid on

p.LineWidth = 2

xlim([-127 127]);
xlabel('Atraso (\tau)','FontSize',19.5,'FontWeight','bold')
ylabel('Autocorrelação R_C (\tau)','FontSize',19.5,'FontWeight','bold')


clear; clc;
x=lfsr([1 1 1 1 0 0 0 1],[1 0 0 0 0 0 0]) %g_1(x)=xˆ5+xˆ4+xˆ2+x+1
y=lfsr([1 0 0 0 0 0 1 1],[1 0 0 0 0 0 0]) %g_2(x)=xˆ5+xˆ4+xˆ3+x+1
N=127; %Period of the sequence N=2ˆL-1
Rxy=sequence_correlation(x,y,-127,127)%cross-correlation

p = plot(-127:1:127,Rxy);   %Plot cross-correlation
grid on
p.LineWidth = 2
xlim([-127 127]);
ylim([-45 45]);
xlabel('Atraso (\tau)','FontSize',19.5,'FontWeight','bold')
ylabel({'Correlação Cruzada R_{12} (\tau)'},'FontSize',19.5,'FontWeight','bold')


clear; clc;
N=127;%period of Gold code
G1=[1 1 1 1 0 0 0 1]; G2=[1 0 0 1 0 0 0 1]; %feedback connections
X1=[1 0 0 0 0 0 0];   X2=[1 0 0 0 0 0 0]; %initial states of LFSRs
y1=gold_code_generator(G1,G2,X1,X2);%Generate Gold code

G1=[1 1 1 1 0 0 0 1]; G2=[1 0 0 1 0 0 0 1]; %feedback connections
X1=[0 0 0 0 0 1 1];   X2=[1 0 0 0 0 0 0]; %initial states of LFSRs
y2=gold_code_generator(G1,G2,X1,X2);%Generate Gold code

Ry1y2= sequence_correlation(y1,y2,-130,130);%cross-correlation
Ry1y1= sequence_correlation(y1,y1,-130,130);%cross-correlation
Ry2y2= sequence_correlation(y2,y2,-130,130);%cross-correlation

p = plot(-130:1:130,Ry1y2);%plot auto-correlation
grid on
p.LineWidth = 2
xlim([-127 127]);
ylim([-20 20]);
xlabel('Atraso (\tau)','FontSize',19.5,'FontWeight','bold')
ylabel({'Correlação Cruzada R_{12} (\tau)'},'FontSize',19.5,'FontWeight','bold')

p = plot(-130:1:130,Ry1y1);%plot auto-correlation
grid on
p.LineWidth = 2
xlim([-127 127]);
%ylim([-20 20]);
xlabel('Atraso (\tau)','FontSize',19.5,'FontWeight','bold')
ylabel({'Autocorrelação R_C (\tau)'},'FontSize',19.5,'FontWeight','bold')

p = plot(-130:1:130,Ry2y2);%plot auto-correlation
grid on
p.LineWidth = 2
xlim([-127 127]);
%ylim([-20 20]);
xlabel('Atraso (\tau)','FontSize',19.5,'FontWeight','bold')
ylabel({'Autocorrelação R_C (\tau)'},'FontSize',19.5,'FontWeight','bold')
