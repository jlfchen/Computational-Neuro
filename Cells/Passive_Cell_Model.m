% Modeling the passive isopotential cell
% Date uploaded: 5/21/20

% (electrical/biological) cellular properties 
r = 10; sa = (r*0.0001).^2*4*pi; % radius [µm], surface area [cm^2]
disp("Surface Area: "+sa +" cm2") % assumption: spherical

Cm = 1; % approximation of capacitance of biological membranes [µF/cm2]
C = Cm*sa; % capacitance [µF]
disp("Capacitance: "+ C +" ÂµF") 

sleak = 0.3; % specific leak conductance [mS/cm2]
gleak = sleak*sa; % membrane conductance [mS]
disp("Membrane conductance: "+gleak)
% corresponding resistance(MOhms)
res = 1/gleak*10.^-3; % scaled to MOhms (originally 1/mS = 1000 Ohms)
disp("Resistance: "+ res +" MOhms")

tao = C./gleak; % membrane time constant (ms)
disp("Tao: "+tao+" ms")

I0 = 10; % injected current [pA]
t0 = 2; t1 = 22; % duration of injection [ms]
ti = 0; tf = 40; % total time frame
Eleak = -68; % resting potential of cell [mV]

syms t
assume (t>ti & t<tf) % Note: I0 pA scaled to muA
% exact solution to passive cell differential equation 
V = Eleak + I0*10.^(-6).*tao./C.*piecewise(...
    t < t0, 0,...
    t0 < t & t < t1, 1-exp(-(t-t0)/tao),...
    t > t1, exp(-(t-t1)./tao)-exp(-(t-t0)./tao));
% Analytical Solution
figure, fplot(V,'Linewidth',1) 

% set time steps and initial conditions
dt = 0.01; t = ti:dt:tf; n = length(t); % [ms]
Vfe = zeros(1,n); Iel = zeros(1,n); 
Vfe(1) = Eleak; Iel(t0/dt+1:t1/dt+1) = I0*10.^(-6); % [muA]
% Forward Euler Approximation
for i = 1:n-1
    Vfe(i+1) = Vfe(i)+(-gleak.*(Vfe(i)-Eleak)+Iel(i)).*dt./C;
end
hold on
% Numerical Solution
plot(t,Vfe, '--','Linewidth',1)
axis([ti tf Eleak -65]), xlabel('t (ms)'), ylabel('V (mV)') 
legend({'V (analytic)','V (numeric)'})
 
% Capacitance and Resistive Currents
Ileak = 10^6.*gleak.*(V-Eleak); % leak (resistive) current
figure, fplot(Ileak,'Linewidth',1)

syms t
assume (t>ti & t<tf)
I0f = piecewise(t<t0, 0, t0<t & t<t1, I0, t>t1,0); % external current
IC = I0f - Ileak; % capacitive current
hold on, fplot(IC,'Linewidth',1)
xlim auto
axis([ti tf -10 Inf]), xlabel('t (ms)'), ylabel('I (pA)')
legend({'Ileak (resistive)','IC (capacitive)'})

!git add Passive_Cell_Model.m
!git commit -m "Model of isopotential cell."
!git push
