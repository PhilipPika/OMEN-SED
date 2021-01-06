%% Supplementary code to replicate figure 1b: Approximation of the RCM with 
% a multi-F approach

% CITATION:Pika et al. (2020), GMD. DOI:
% Written for Matlab 2020a by Philip Pika, ver. 2018-20-10
% Contact: philip.pika@ulb.be

%% Begin Code
clearvars; close all; clc
set(0,'defaultLineLineWidth', 2)
set(0,'DefaultAxesFontSize',18)
figure
% define parameters
nG=400;                 % number of G fractions used to approximate a-nu continuum
nu=0.125;%0.125;               % nu
a=3e-4;%3e-4;                    % a
POC0=2;                 % TOC at SWI (%)
w0=0.074;%0.004;        % sed rate in cm/yr
db0=27.6;               % bioturbation coefficient in cm2/yr
zinf = 100;
dz=10;                  % resolution
emin= log10(10e-16);%...
    %nu./(a+(zinf./w0))...
    %) - 1;                    % lower k limit for k-bins of multi-G approximation, i.e. k=1e-15 yr-1
% subtracted 1 to be on the conservative side.
emax=-log10(a)+2;       % upper k limit for k-bins of multi-G approximation
if emin >= emax;error('emin >= emax, this cannot be!');end

%__________________________________________________________________________
%   multi-G approximation of the reactive continuum model for the
%   bioturbated zone
%   (analytical solution of multi-G advection-diffusion-reaction equation,
%   k and F of the multi-G model are calculated by integrating the reactive
%   continuum model)
%__________________________________________________________________________

% calculate first and last bin of the multi-G model from a-nu distribution
k(1)= 10^(emin);
kk(1)=10^(emin);
F(1) = gammainc(a*10^emin,nu,'lower');
kk(nG)=10^(emax);
k(nG)=10^(emax);
F(nG) = gammainc(a*10^emax,nu,'upper');
% calculate intemediate bins of the multi-G model from a-nu distribution
G=2:nG-1;
ne=emin+(1:nG-2).*(emax-emin)./(nG-1);
kk(2:nG-1)=10.^ne; clear ne
G_inc_0 = gammainc(a*kk(1:nG-2),nu,'upper');
G_inc_1 = gammainc(a*kk(2:nG-1),nu,'upper'); % G = 2:end-1
F(G) = (G_inc_0 - G_inc_1);
k(G)=kk(1:nG-2)+(kk(2:nG-1)-kk(1:nG-2))/2; clear G*
if abs(sum(F)-1) > 0.0001;warning('F~=1!!');end

Fnonbioi = F.* ( POC0*(1-0.7)*w0 );
POC0i = F.*POC0;
sum(POC0i)

a1 = (w0-(w0^2+4*db0*k).^(1/2))/(2*db0);
b1 = (w0+(w0^2+4*db0.*k).^(1/2))/(2*db0);
a2=(-k./w0);

kapp_100G = nu/a


plot(log10(k), POC0i,'b')
hold on
bar(log10(k), POC0i)
xlabel('log_{10} (OM k) [yrs^{-1}]')
ylabel('Fraction of total OM')
box on;
%title(['Sum all F: '  num2str(sum(F),20)])
%ntitle(['emin: ' num2str(emin) ', emax: ' num2str(emax)])
set(gcf,'Position',[680 678 860 420],'color','w')
export_fig(gcf,'ApproximationOfRCM','-pdf','-png','-r150')
