% Physics567 final presentation fig7 run for thick blue curve
%dendritic injection only
% time vector in ms

%% Parameters
RT = 65; % in Mohm
RS = 50; % in Mohm
RD = 43; % in Mohm
CS = 0.26; % in nF
CD = 0.12; % in nF
Vrs = -70; % in mV
Vrd = -60; % in mV
theta = -47; %in mV (soma threshold)
reset = -52; % in mV (soma reset)

% potassium current
gAHP = 4E-3 ; % in mS
EK = -90 ;% in mV
tauK = 80 ;% in ms

%calcium current
gCa = 70E-3 ;% in mS
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
tau = 3; %ms


%test soma threshold
mus = 0.25:0.05:0.8;
%injd
%mus = 1.0:0.03:1.33;
            
%% iterations
cvs=[];
num_iters = 5;
freqs_sum = zeros(1, length(mus));
for i = 1:num_iters
    j=0;
    %generate noisy current imput
    currents1 = [];
    muvec = [];
    for mu = mus
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus))
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
    for mu = mus
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus))
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
    %Iinjd = currents2;
    Iinjd = zeros(1, length(currents2));
    
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
    window_size = floor(length(tvec) / length(mus)); %step duration of 2s
    num_spikes_s = []; %soma
    cv = [];
    loc_s = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:length(mus)-1
        num_spikes_s(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        loc_s{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc_s{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1});%getCV(isi{i+1});
    end
    freqs_s = ((num_spikes_s * 1000) ./ (window_size));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs_s;
end
ave_freqs_tb = freqs_sum ./ num_iters; %soma 

%store data for thin blue curve
vs_tb = vs; %thin blue
vd_tb = vd;
ICa_tb = ICa;
Iahp_tb = Iahp;
mus_tb =mus;