%% Supplementary code to replicate figure 7: Dependency of SMTZ with 
% parameter a and sedimentation rate

% CITATION:Pika et al. (2020), GMD. DOI:
% Written for Matlab 2020a by Philip Pika, ver. 2018-20-10
% Contact: philip.pika@ulb.be

%% Begin Code
close all;clc; clearvars
addpath ../OMEN_MATLAB
% -------------------------------------------------------------------------
%% START SIMULATIONS HERE
% -------------------------------------------------------------------------

% Defining a and nu
A =  linspace(-1,8,200);
NU = 0.2;
anu = combvec(A , NU)';

SR = linspace(benthic_main.sedrate(50), benthic_main.sedrate(250),10)./1; % depth in [m]
hC = linspecer(size(SR,2)); % requires linspcer.m see mathworks.com
OutputShelf = ones(size(anu,1),10)*NaN;
tic
for idx_SR = 1:size(SR,2)
    for idx_anu = 1:size(anu,1)
        BCin = [...
            200 ...         % nanmedian(abs(double(Shelf.Depth)))
            270.55 ...      % nanmedian(Shelf.Oxygen)
            8.29  ...       % nanmedian(Shelf.Nitrate)
            .8  ...
            1.46 ...        % nanmedian(Shelf.Temp)
            843.89e-003 ... % nanmedian(Shelf.PO4)
            500 ...         % nG
            SR(idx_SR)];
        if isnan(BCin) == 1;break;end
        
        [OMEN_Output] = initiateOMEN_SED_RCM(BCin,anu(idx_anu,:));
        
        OutputShelf(idx_anu,:) = [...
            OMEN_Output.flxO2D_PP ...
            OMEN_Output.flxO2D_PP*3.65e5...
            OMEN_Output.zox ...
            OMEN_Output.flxswiNO3...
            OMEN_Output.zno3...
            OMEN_Output.flxswiNH4...
            OMEN_Output.flxswiSO4 ...
            OMEN_Output.zso4...
            OMEN_Output.swi.C0 / (1e-2/12*OMEN_Output.bsd.rho_sed) ...
            OMEN_Output.flxswiO2];
    end
    Output{idx_SR} = OutputShelf;
end
toc

%%
hold on
for idx_SR = 1:size(SR,2)
    plot(A,-Output{:,idx_SR}(:,8)./100,'Color',hC(idx_SR,:))
end

hold off
xlabel('log_{10}(Parameter a) [yrs]')
ylabel('SMTZ [m]')
box on; set(gcf, 'Color', 'w'); lighting phong;
hleg = legend(strsplit(num2str(round(SR,3))),'Location','best');
title(hleg, 'Sedimentation rates [cm/yr]'); set(gcf,'Position', [1160 991 1073 527])
xlim([0 7])
export_fig(gcf,['parameter_a_vs_SMTZ_TOC=' num2str(OMEN_Output.swi.C0*12*100/2.5)], '-png','-m3.4')
%%
function [OMEN_Output] = initiateOMEN_SED_RCM(BCin,anu)
sD.SFD = BCin(1);
sD.dissO2 = BCin(2);
sD.NO3 = BCin(3);
sD.TOC = BCin(4);
sD.Temp = BCin(5);
sD.PO4 = BCin(6);
sD.SR = BCin(8); % OR NaN
sD.nG = BCin(7);

n = 1;

OMEN_Output = benthic_test.OMEN_with_DOU_input(sD,anu(n,:));
end