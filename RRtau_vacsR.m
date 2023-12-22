
clear all;
clc;

hh = 0.05;
TAUEND = 15; 
Rtaui = 21/hh;

TAU = 28:3.5:7*TAUEND; % waiting period for the 2nd dose: 1 weeks to 12 weeks
EPS1 = 0.3:0.025:0.6; % efficacy of the 1st dose against infection
EPS2 = 0.75:0.15:0.9; % efficacy of the 2nd dose against infection, 0.75 and 0.9

m = length(TAU);
e1 = length(EPS1);
e2 = length(EPS2);

RdeathP = zeros(e1, e2, m);
RdeathR = zeros(e1, e2, m);
RhospP = zeros(e1, e2, m);
RhospR = zeros(e1, e2, m);
RinciP = zeros(e1, e2, m);
RinciR = zeros(e1, e2, m);


DV = 400; % 400 % maxinum number of doses everyday, D
tend = 105; % 100; % 105 for RR during the delayed tau days, 100 for the first 100 days
eend = tend/hh;
     for ee1 = 1:e1
         epsilon1 = EPS1(ee1);
         for ee2 = 1:e2
             epsilon2 = EPS2(ee2);
            [Pinci, Sinci, inciS, inciV1, inciV1w, inciE1, inciE1w, inciV2, inciV2w, inciW2, inciW2w, inciR1w, inciR2w, Dp, Ds, Phosp, Shosp, nu1, nu2, nu2w] = vacsR(hh, DV, epsilon1, epsilon2, tend);
            for ii = 1:m
                 tau = TAU(ii);
                 taui = tau/hh;
              
                 RdeathP(ee1, ee2, ii) = Dp(taui)-Dp(Rtaui);
                 RdeathR(ee1, ee2, ii) = Ds(taui)-Ds(Rtaui);
                 RhospP(ee1, ee2, ii) = sum(Phosp(Rtaui:taui))*hh;
                 RhospR(ee1, ee2, ii) = sum(Shosp(Rtaui:taui))*hh;
                 RinciP(ee1, ee2, ii) = sum(Pinci(Rtaui:taui))*hh;
                 RinciR(ee1, ee2, ii) = sum(Sinci(Rtaui:taui))*hh;
             end
         end
     end

 save('R1_400.mat', "RinciR","RinciP", "RhospR","RhospP","RdeathR","RdeathP")
%save('RR100_R1400.mat', "RinciR","RinciP", "RhospR","RhospP","RdeathR","RdeathP")


%%
% Simulation for optimal dosing interval with waning parameter set I/II

function [Pinci, Sinci, inciS, inciV1, inciV1w, inciE1, inciE1w, inciV2, inciV2w, inciW2, inciW2w, inciR1w, inciR2w, Dp, Ds, Phosp, Shosp, nu1, nu2, nu2w] = vacsR(hh, DV, epsilon1, epsilon2, tend)

Rtau = 21; % recommanded dosing interval after first dose

% parameters with waning immunity set II:

% omega1 = 1/90; % duration of full protection following recovery (R)
% omega2 = 1/60; % duration of partial protecttion following recovery (R1w)
% xi1 =  1/60; % duration of protection after first dose
% xi2 =  1/90; % duration of protection after second dose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters with waning immunity set I:

omega1 = 1/120; % duration of full protection following recovery (R)
omega2 = 1/90; % duration of partial protecttion following recovery (R1w)
xi1 =  1/90; % duration of protection after first dose
xi2 =  1/120; % duration of protection after second dose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NP = 100000; % initial population size
%DV = 500; % maximum daily doses;
kk = 0.9; % proportion of D for second dose
TV = 0.8*NP; % total number of vaccine

%--fixed parameters

gamma = 1/5.5; % infectious period 
mu = 1/10; % average length of hospital stay
d = 0.125; % proportion of hospitalised individuals who die

%epsilon1 = 0.3; % effectiveness of 1st dose against infection 30% - 60%
%epsilon2 = 0.7; % effectiveness of 2st dose against infection 60% - 90%

h = 1/3;  % duration of infectious period before hospitalisation 2-4 days
sigma = 0.0138; % chosed such that alpha = 0.025
sigma1 = 0.00342; % chosen such that delta1 = 0.75
sigma2 = 0.000682; % chosen such that delta2 = 0.95

%alpha = 0.025; % fraction of hospitalized infected individuals
%delta1 = 0.75; % 1st dose agaist hospitalization 60% - 90%
%delta2 = 0.95; % 2nd dose against hospitalization fixed

eta1 = 0.75; % effectiveness of 1st dose against death
eta2 = 0.95; % effectiveness of 2nd dose against death

R0 = 1.10; % 1.8; % reprodution number
beta=R0*gamma/NP; % transmission rate

%--set up simulation 
%hh = 0.05; % step size for simulation
T0 = 0;  % start time
TN = tend;  % end time % 100 or 105
tt = T0:hh:TN;
N = length(tt);  % number of steps

Rtaui = Rtau/hh;

%-----------initialize variables
S = zeros(1,N); % susceptible
V1 = zeros(1, N);
V1w = zeros(1, N);
E1 = zeros(1, N);
E1w = zeros(1, N);
V2 = zeros(1, N);
V2w = zeros(1, N);
W2 = zeros(1, N);
W2w = zeros(1, N);
I = zeros(1,N);
I1p = zeros(1,N);
I1s = zeros(1, N);
I2p = zeros(1,N);
I2s = zeros(1, N);
H = zeros(1,N);
H1p = zeros(1,N);
H1s = zeros(1, N);
H2p = zeros(1,N);
H2s = zeros(1, N);
Dp = zeros(1,N);
Ds = zeros(1, N);
R = zeros(1,N);
R1w = zeros(1, N);
R2w = zeros(1, N);

vttau = zeros(1, N);
vwttau = zeros(1, N);
Lambda = zeros(1, N);


cumnu = zeros(1, N);
nu = zeros(1,N);
nu1 = zeros(1,N); 
Tnu2 = zeros(1,N);
nu2 = zeros(1, N);
nu2w = zeros(1, N);

Pinci = zeros(1, N);
Sinci = zeros(1, N);

inciS = zeros(1,N);
inciV1 = zeros(1,N);
inciV1w = zeros(1, N);
inciE1 = zeros(1, N);
inciE1w = zeros(1, N);
inciV2 = zeros(1,N);
inciV2w = zeros(1, N);
inciW2 = zeros(1, N);
inciW2w = zeros(1, N);
inciR1w = zeros(1, N);
inciR2w = zeros(1, N);


Phosp = zeros(1,N);
Shosp = zeros(1, N);

% initial conditions
I(1) = 10;
S(1) = 95000 - I(1);
R(1) = 3000;
R1w(1) = 2000;

%--main loop--Non-standard method
for i = 1 : (N-1)
    %--calculate incidence
    Lambda(i) = I(i) + I1p(i) + I1s(i) + I2p(i) + I2s(i);
    inciS(i) = beta *S(i)*Lambda(i);
    inciV1(i) = beta *(1-epsilon1)*V1(i)*Lambda(i);
    inciV1w(i) = beta*V1w(i)*Lambda(i);
    inciE1(i) = beta*(1-epsilon1)*E1(i)*Lambda(i);
    inciE1w(i) = beta*E1w(i)*Lambda(i);
    inciV2(i) = beta *(1-epsilon2)*V2(i)*Lambda(i);
    inciV2w(i) = beta*(1-epsilon1)*V2w(i)*Lambda(i);
    inciW2(i) = beta*(1-epsilon1)*W2(i)*Lambda(i);
    inciW2w(i) = beta*W2w(i)*Lambda(i);
    inciR1w(i) = beta*(1-epsilon1)*R1w(i)*Lambda(i);
    inciR2w(i) = beta*R2w(i)*Lambda(i);

    Pinci(i) = inciS(i)+inciV1(i)+inciV1w(i)+inciE1(i)+inciE1w(i)+inciV2(i)+inciV2w(i)+inciW2(i)+inciW2w(i); % incidence of primary infection
    Sinci(i) = inciR1w(i)+inciR2w(i); % incidence of secondary infection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if i > Rtaui
        vwttauh = zeros(1, Rtaui+1);
     for uu = 1:Rtaui
          intLambda1(i-Rtaui) = 0;
          intLambda2(i-uu) = 0;
        for ss = (i-Rtaui):(i)
            intLambda1(ss+1) = intLambda1(ss) + Lambda(ss)*hh;
        end
        vttau(i) = nu1(i-Rtaui)*exp(-xi1*Rtau - beta*(1-epsilon1)*intLambda1(i+1));
        for ss2 = (i-uu):i
            intLambda2(ss2+1) = intLambda2(ss2) + Lambda(ss2)*hh;
        end
        vwttauh(uu+1) = vwttauh(uu) +xi1*nu1(i-Rtaui)*exp(-xi1*(Rtaui-uu)*hh)*exp(-beta*(1-epsilon1)*(intLambda1(i)-intLambda2(i-uu+1)))*exp(-beta*intLambda2(i-uu+1))*hh;
     end 
     vwttau(i) = vwttauh(Rtaui+1);
    end
    
    % define nu2
    if i > Rtaui+1
        if S(i)>0
            Tnu2(i) = min(kk*DV, E1(i)+E1w(i));
            nu1(i) = DV - Tnu2(i);
            nu2(i) = Tnu2(i)*E1(i)/(E1(i) + E1w(i));
            nu2w(i) = Tnu2(i)*E1w(i)/(E1(i)+ E1w(i));
        else
            Tnu2(i) = min(DV, E1(i)+E1w(i));
            nu1(i) = 0;
            nu2(i) = Tnu2(i)*E1(i)/(E1(i) + E1w(i));
            nu2w(i) = Tnu2(i)*E1w(i)/(E1(i)+ E1w(i));
        end
    end

    if i < Rtaui + hh % in the first 0 - tau days, no one receive the 2nd dose
        nu1(i) = DV;
        nu2(i) = 0;
        nu2w(i) = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cumnu(i) > TV
        nu1(i) = 0;
        nu2(i) = 0;
        nu2w(i) = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
    S(i+1) = (S(i) - hh*nu1(i))/(1+hh*beta*Lambda(i));
    V1(i+1) = (V1(i)+ hh*nu1(i)-hh*vttau(i))/(1+hh*beta*(1-epsilon1)*Lambda(i)+hh*xi1);
    V1w(i+1) = (V1w(i) + hh*xi1*V1(i) - hh*vwttau(i))/(1+hh*beta*Lambda(i));
    E1(i+1) = (E1(i) + hh*vttau(i)-hh*nu2(i))/(1+hh*beta*(1-epsilon1)*Lambda(i)+hh*xi1);
    E1w(i+1) = (E1w(i) + hh*vwttau(i) + hh*xi1*E1(i)-hh*nu2w(i))/(1+hh*beta*Lambda(i));
    V2(i+1) = (V2(i) + hh*nu2(i))/(1+ hh*beta*(1-epsilon2)*Lambda(i) + hh*xi2);
    V2w(i+1) = (V2w(i) + hh*nu2w(i))/(1+hh*beta*(1-epsilon1)*Lambda(i)+hh*xi1);
    W2(i+1) = (W2(i) + hh*xi2*V2(i))/(1+hh*beta*(1-epsilon1)*Lambda(i));
    W2w(i+1) = (W2w(i) + hh*xi1*V2w(i))/(1+hh*beta*Lambda(i));    
    I(i+1) = (I(i) + hh*inciS(i))/(1+hh*((1-sigma)*gamma+sigma*h));
    I1p(i+1) = (I1p(i) + hh*(inciV1(i)+inciE1(i) + inciV1w(i)+ inciE1w(i) + inciW2w(i)))/(1+hh*((1-sigma1)*gamma+sigma1*h)); % seperate the primary/secondary infection
    I2p(i+1) = (I2p(i) + hh*(inciV2(i)+ inciV2w(i) + inciW2(i)))/(1+hh*((1-sigma2)*gamma+sigma2*h));
    I1s(i+1) = (I1s(i) + hh*inciR2w(i))/(1+hh*((1-sigma1)*gamma+sigma1*h));
    I2s(i+1) = (I2s(i) + hh*inciR1w(i))/(1+hh*((1-sigma2)*gamma+sigma2*h));

    H(i+1) = (H(i)+hh*sigma*h*I(i))/(1+hh*mu);
    H1p(i+1) = (H1p(i)+hh*sigma1*h*I1p(i))/(1+hh*mu); % seperate hospitalization from primary/secondary infection
    H2p(i+1) = (H2p(i)+hh*sigma2*h*I2p(i))/(1+hh*mu);
    H1s(i+1) = (H1s(i)+hh*sigma1*h*I1s(i))/(1+hh*mu);
    H2s(i+1) = (H2s(i)+hh*sigma2*h*I2s(i))/(1+hh*mu);

    Dp(i+1) = Dp(i) + hh*(d*mu*H(i) + d*mu*(1-eta1)*H1p(i) + d*mu*(1-eta2)*H2p(i));
    Ds(i+1) = Ds(i) + hh*(d*mu*(1-eta1)*H1s(i) + d*mu*(1-eta2)*H2s(i));
    R(i+1) = (R(i) + hh*((1-sigma)*gamma*I(i) + (1-sigma1)*gamma*(I1p(i)+I1s(i)) + (1-sigma2)*gamma*(I2p(i)+I2s(i))...
        + (1-d)*mu*H(i) + (1-(1-eta1)*d)*mu*(H1p(i)+H1s(i)) +(1-(1-eta2)*d)*mu*(H2p(i)+H2s(i))) )/(1+hh*omega1);
    R1w(i+1) = (R1w(i) + hh*omega1*R(i) )/(1+ hh*omega2 + hh*beta*(1-epsilon1)*Lambda(i));
    R2w(i+1) = (R2w(i) + hh*omega2*R1w(i))/(1+hh*beta*Lambda(i));

    nu(i) = nu1(i)+ nu2(i)+nu2w(i);
   
    cumnu(i+1) = cumnu(i) + nu(i)*hh;
    Phosp(i) = sigma*h*I(i) + sigma1*h*I1p(i) + sigma2*h*I2p(i); % incidence of hosp from primary infection
    Shosp(i) = sigma1*h*I1s(i) + sigma2*h*I2s(i);
end

end


