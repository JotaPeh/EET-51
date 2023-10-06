clear; clc;
y=lfsr([1 0 0 0 0 0 1 1],[0 0 0 0 0 0 1]) %g(x)=xˆ5+xˆ2+1
N=127;%Period of the sequence N=2ˆL-1
Ryy=sequence_correlation(y,y,-130,130)%normalized autocorr.
plot(-130:1:130,Ryy) %Plot autocorrelation
grid on
xlim([-130 130]);

x=lfsr([1 1 1 1 0 0 0 1],[0 0 0 0 0 0 1]) %g_1(x)=xˆ5+xˆ4+xˆ2+x+1
y=lfsr([1 0 0 0 0 0 1 1],[0 0 0 0 0 0 1]) %g_2(x)=xˆ5+xˆ4+xˆ3+x+1
N=127; %Period of the sequence N=2ˆL-1
Rxy=sequence_correlation(x,y,-127,127)%cross-correlation
plot(-127:1:127,Rxy)%Plot cross-correlation
grid on
xlim([-127 127]);

clear; clc;
G1=[1 1 1 1 0 0 0 1]; G2=[1 0 0 1 0 0 0 1]; %feedback connections
X1=[0 0 0 0 0 0 1];   X2=[0 0 0 0 0 0 1]; %initial states of LFSRs
y1=lfsr(G1,X1); y2=lfsr(G2,X2);N=127; %m-sequence 1 and 2
Ry1y2= sequence_correlation(y1,y2,0,127);%cross-correlation
plot(0:1:127,Ry1y2)%plot correlation
grid on
xlim([0 127]);

N=127;%period of Gold code
G1=[1 1 1 1 0 0 0 1]; G2=[1 0 0 1 0 0 0 1]; %feedback connections
X1=[0 0 0 0 0 0 1];   X2=[0 0 0 0 0 0 1]; %initial states of LFSRs
y=gold_code_generator(G1,G2,X1,X2);%Generate Gold code
Ryy= 1/N*sequence_correlation(y,y,0,127);%auto-correlation
plot(0:1:127,Ryy);%plot auto-correlation
grid on
xlim([0 127]);
