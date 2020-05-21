% Effect of A-type (potassium) current on spike latency 
% model of stellate cells of cerebellum (type of inhibitory neuron)
% Date uploaded: 5/21/20

% steady-state activation and inactivation functions of channels
% A (potassium) channel
ass = @(V) (1+exp(-(V+27)/8.8)).^(-1); % activation
bss = @(V) (1+exp((V+68)/6.6)).^(-1); % inactivation
V = linspace(-100,20);
figure, plot(V,ass(V)), hold on, plot(V,bss(V))
xlabel('V (mV)'), text(-70,0.7,'binf'), text(-20,0.5,'ainf')
% Na channel
mss = @(V) (1+exp(-(V+35)/4)).^(-1); % activation
hss = @(V) (1+exp((V+35)/4)).^(-1); % inactivation
figure, plot(V,mss(V)), hold on, plot(V,hss(V))
xlabel('V (mV)'), text(-45,0.65,'hinf'), text(-20,0.9,'minf')

% inactivation time constant of Na channel
taoh = @(V) 12992./(4*pi*(V+74).^2+784)-0.15;
figure, plot(V,taoh(V))
ylabel('taoh (ms)'), xlabel('V (mV)') %'\tao_{h} (ms)'

% cellular constants (biological parameters)
taob = 15; taon = 0.5; nss = mss; % time constant [ms]
gNa = 30; gK = 7; gleak = 1; gA = 16; % [mS/cm^2]
VNa = 45; VK = -90; Vleak = -70; %[mV]
C = 1.5; %[µF/cm2]

% model parameters
I0 = 24; t0 = 5; t1 = 7; % current injection [muA/cm^2] [ms]
ti = 0; tf = 20; dt = 0.005; % time frame [ms]
Vrest = -69.394; % resting potential [mV]
N = (tf-ti)./dt+1; t = ti:dt:tf;
% initializations
Vfe = zeros(1,N); IK = zeros(1,N-1); INa = zeros(1,N-1);
n = zeros(1,N); m = zeros(1,N); h = zeros(1,N);
Iel = zeros(1,N); Iel(t0/dt+1:t1/dt+1) = I0; % [muA] *10.^(-6)
Vfe(1) = Vrest;
n(1) = nss(Vfe(1)); m(1) = mss(Vfe(1)); h(1) = hss(Vfe(1)); 
% Forward Euler (single compartment model)
for i = 1:N-1
    V = Vfe(i);
    n(i+1) = n(i)-((n(i)-mss(V))/taon)*dt;
    h(i+1) = h(i)-((h(i)-hss(V))/taoh(V))*dt;
    IK(i) = gK.*n(i+1)*(V-VK);
    INa(i) = gNa.*mss(V).*h(i+1)*(V-VNa);    
    Vfe(i+1) = V-(gleak.*(V-Vleak)+IK(i)+INa(i)-Iel(i)).*dt./C;
end
figure, plot(t,Vfe)

% above model with A current
Vrest = -70.837; % [mV]
Vafe = zeros(1,N); IaK = zeros(1,N-1); IaNa = zeros(1,N-1); Ia = zeros(1,N-1);
Vafe(1) = Vrest; 
na = zeros(1,N); ma = zeros(1,N); ha = zeros(1,N); b = zeros(1,N);
na(1) = nss(Vafe(1)); ma(1) = mss(Vafe(1)); ha(1) = hss(Vafe(1)); b(1) = bss(Vafe(1));
for i = 1:N-1
    V = Vafe(i);
    na(i+1) = na(i)-((na(i)-mss(V))/taon)*dt;
    ha(i+1) = ha(i)-((ha(i)-hss(V))/taoh(V))*dt;
    b(i+1) = b(i)-((b(i)-bss(V))/taob)*dt;
    IaK(i) = gK.*na(i+1)*(V-VK);
    IaNa(i) = gNa.*mss(V).*ha(i+1)*(V-VNa);
    Ia(i) = gA*ass(V)*b(i+1)*(V-VK);
    Vafe(i+1) = V-(Ia(i)+gleak.*(V-Vleak)+IaK(i)+IaNa(i)-Iel(i)).*dt./C;
end
hold on, plot(t,Vafe), xlim([0,15])
ylabel('V (mV)'), xlabel('t (ms)'), legend('without A','with A')

% current flowing through various channels
% without A current
figure, plot(t(1:end-1),IK)
hold on, plot(t(1:end-1),INa)
hold on, plot(t(1:end-1),zeros(1,N-1))
xlim([6,9]), ylabel('muA/cm2'), xlabel('t (ms)'), legend('IK','INa','IA')
% with A current
figure, plot(t(1:end-1),IaK)
hold on, plot(t(1:end-1),IaNa)
hold on, plot(t(1:end-1),Ia)
xlim([6,9]), ylabel('muA/cm2'), xlabel('t (ms)'), legend('IK','INa','IA')

!git add A_Current.m
!git commit -m "Model of stellate cells of cerebellum."
!git push