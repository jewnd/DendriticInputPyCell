% Physics567 final presentation fig7 run for thin blue curve

%% run for soma + 750pA dendrite daya
clear;clc;
global RT RS RD CS CD Vrd Vrs theta Eca reset

%% Parameters
RT = 65; % in Mohm %model distal dendrite (original 65)
RS = 50; % in Mohm
RD = 43; % in Mohm
CS = 0.26; % in nF
CD = 0.12; % in nF
Vrs = -70; % in mV
Vrd = -60; % in mV
theta = -47; %in mV (soma threshold)
reset = -52; % in mV (soma reset)

% potassium current
gAHP = 4E-3 ; % in uS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current (in nA)
gCa = 70E-3 ;% in uS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 24000; 
dt = 0.1;
tvec = t1:dt:t2;

%vectors (vs,vd,m,h,Iahp,ICa)
vs = zeros(size(tvec));
vd = zeros(size(tvec));
vs(1) = Vrs;
vd(1) = Vrd;
m = ones(size(tvec));
h = zeros(size(tvec));
Iahp = zeros(size(tvec));
ICa = zeros(size(tvec));

Eca = 137;  % Equilibrium potential for Ca 
            % [https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html]

%% Noisy input current

sigma = 0.3; %nA
%sigma = 0;
tau = 3; %ms

% mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
%     0.18, 0.21, 0.24, 0.27, 0.30, 0.33]; %in nA
% test soma threshold
d_inj = 0.75; %in nA
mus_s = 0:0.03:0.33;%step current somatic injection
mus_d = ones(1,length(mus_s)).*d_inj; %steady dendritic injection
% injd
% mus = 1:0.06:2.3;
            
%% iterations
cvs=[];
num_iters = 10;
freqs_sum = zeros(1, length(mus_s));
for i = 1:num_iters
    j=0;
    %generate noisy current imput
    currents1 = [];
    muvec = [];
    for mu = mus_s
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus_s))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
           
        end
        muv = ones(size(current)).*mu;
        currents1 = horzcat(currents1, current);
        muvec = horzcat(muvec,muv);
    end
    
    for i = length(currents1):length(tvec)
        currents1(i) = mu(end);
    end
    movave1 = movmean(currents1, 1000);

    currents2 = [];
    for mu = mus_d
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus_d))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
        end
        currents2 = horzcat(currents2, current);
    end
    for i = length(currents2):length(tvec)
        currents2(i) = mu(end);
    end
    movave2 = movmean(currents2, 1000);
    
    Iinjs = currents1;
    Iinjd = currents2;
    %Iinjd = zeros(1, length(currents1));
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));  
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            vs(i+1) = reset;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            vs(i+1) = 10; %for 1ms
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus_s)); %step duration of 2s
    num_spikes_s = []; %soma
    cv = [];
    loc_s = cell(size(mus_s));
    isi = cell(size(mus_s));
    for i = 0:length(mus_s)-1
        num_spikes_s(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        loc_s{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc_s{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1});%getCV(isi{i+1});
    end
    freqs_s = ((num_spikes_s * 1000) ./ (window_size));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs_s;
end
ave_freqs_sb = freqs_sum ./ num_iters; %soma 

%store data for thin blue curve
vs_sb = vs; %thin blue
vd_sb = vd;
ICa_sb = ICa;
Iahp_sb = Iahp;
mus_sb =mus_s;


%% plotting

%fig 7B
figure(2)
plot(mus_sb.*1000,ave_freqs_sb,'b-','Linewidth',1.5)
hold on
plot(mus_tb.*1000,ave_freqs_tb,'b-','Linewidth',2.5)
plot(mus_tr.*1000,ave_freqs_tr,'r-','Linewidth',2.5)
plot(mus_sr.*1000,ave_freqs_sr,'r-','Linewidth',1.5)
hold off
xlabel('Somatic and dendritic current means (pA)')
ylabel('Mean firing rates (AP/s)')
text(1300,19, 'n = 10', 'HorizontalAlignment','center')
legend('s+d','s','d','d-Ca2+','Location','Northwest')
xlim([-500 1500])
xticks([-500 0 500 1000 1500])
title('Mean firing rates for somatic and dendritic current injections')

%fig7C1 - 1s
figure(3)
plot(tvec((11*window_size):(11*window_size+10000)), vs_tb((11*window_size):(11*window_size+10000) ),'k-')
hold on
plot(tvec((11*window_size):(11*window_size+10000)), vd_tb((11*window_size):(11*window_size+10000) ),'r-')
plot([22100 22100],[-60 -70],'-k',[22150 22250],[-65 -65],'-k','Linewidth',2)
text(22050,-65,'10 mV','HorizontalAlignment','center')
text(22200, -62,'100 ms','HorizontalAlignment','center')
set(gca,'Visible','off')
hold off
legend('Somatic voltage','Dendritic Voltage','Location','Southeast')

figure(4)
plot(tvec((11*window_size):(11*window_size+10000)), ICa_tb((11*window_size):(11*window_size+10000) ),'r-','Linewidth',1.3)
hold on
plot(tvec((11*window_size):(11*window_size+10000)), Iahp_tb((11*window_size):(11*window_size+10000) ),'b-')
plot([21900 21900],[0 -0.5],'-k',[22150 22250],[-0.6 -0.6],'-k','Linewidth',2)
text(21850,-0.25,'0.5 nA','HorizontalAlignment','center')
text(22200, -0.54,'100 ms','HorizontalAlignment','center')
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) 0.3*tmp(4)])
hold off
set(gca,'Visible','off')
legend('I_C_a','I_A_H_P','Location','Southeast','NumColumns',2)

%fig7C2
figure(5)
plot(tvec((11.4*window_size):(11.4*window_size+10000)), vs_tr((11.4*window_size):(11.4*window_size+10000) ),'k-')
hold on
plot(tvec((11.4*window_size):(11.4*window_size+10000)), vd_tr((11.4*window_size):(11.4*window_size+10000) ),'r-')
plot([22900 22900],[-60 -70],'-k',[22950 23050],[-65 -65],'-k','Linewidth',2)
text(22850, -65,'10 mV','HorizontalAlignment','center')
text(23000,-62,'100 ms','HorizontalAlignment','center')
set(gca,'Visible','off')
hold off
legend('Somatic voltage','Dendritic Voltage','Location','Southeast')

figure(6)
plot(tvec((11.4*window_size):(11.4*window_size+10000)), ICa_tr((11.4*window_size):(11.4*window_size+10000) ),'r-')
hold on
plot(tvec((11.4*window_size):(11.4*window_size+10000)), Iahp_tr((11.4*window_size):(11.4*window_size+10000) ),'b-')
plot([22900 22900],[0.1 1.1],'-k',[22950 23050],[-0.7 -0.7],'-k','Linewidth',2)
text(22870,0.55,'1 nA','HorizontalAlignment','center')
text(23000, -0.6,'100 ms','HorizontalAlignment','center')
% tmp=get(gca,'position');
% set(gca,'position',[tmp(1) tmp(2) tmp(3) 0.3*tmp(4)])
hold off
set(gca,'Visible','off')
legend('I_C_a','I_A_H_P','Location','Southeast','NumColumns',2)


