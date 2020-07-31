% OMEN-SED 1.0 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al (2017) GMD paper

% benthic_test.m
% functions to run OMEN-SED and plot the results

% Command to run the model: benthic_test.run_OMEN
% 1) default sediment-water interface boundary conditions are set as
%    prescribed in default_swi()
% 2) the subroutines for the different properties are called in test_benthic(1,swi):
% 3) results saved in res and are plotted with plot_column()

classdef benthic_test
    % test cases for benthic layer model
    
    properties
    end
    
    methods(Static)
        
        function swi = default_swi()
            % set default SWI conditions
            
            format longEng
            
            bsd = benthic_main();
            %bottom water concentrations
            swi.T = 8.0;                                        % temperature (degree C)
            
            swi.Nitrogen=true;                                  % calculate N (true/false)
            swi.Iron=true;                                      % calculate Fe (true/false)

% WITH TOC CONCENTRATION
            swi.C01_wtprc = 1.63;                                % TOC1 concentration at SWI (wt%)
            swi.C02_wtprc = 1.05;                                % TOC2 concentration at SWI (wt%)
            swi.C01nonbio= swi.C01_wtprc*1e-2/12*bsd.rho_sed;             % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.C02nonbio= swi.C02_wtprc*1e-2/12*bsd.rho_sed;             % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.Fnonbio1 = swi.C01nonbio*(1-bsd.por)*bsd.w;     % Dale: 1.460E-04/2; %    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.Fnonbio2 = swi.C02nonbio*(1-bsd.por)*bsd.w;     % Dale: 1.460E-04/2;  %   % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.C01 = swi.C01nonbio;                            % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.C02 = swi.C02nonbio;                            % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
% WITH TOC FLUX
%            swi.Fnonbio1 = 6.0E-005;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
%            swi.Fnonbio2 = 2.8E-006;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
%            swi.C01nonbio = swi.Fnonbio1/((1-bsd.por)*bsd.w);
%            swi.C02nonbio = swi.Fnonbio2/((1-bsd.por)*bsd.w);
%            swi.C01 = swi.C01nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
%            swi.C02 = swi.C02nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m

%            swi.FeIII0=2.8E-005; %3.0E-006;                   	% FeIII concentration at SWI (mol/cm^3) --> TODO: needs to be a flux!
            Fe_influx = 0.1667*1110.0E-006;                     % FeIII influx from Dale but just 16.67% is available for DIR (See his Tab. 2, footnote d)
            swi.Flux_FeIII0 =  (Fe_influx)*365/100^2;           % Dale 1110 mumol/(m^2 day)   -->  mumol/(cm^2 yr):     *365/100^2
            swi.FeIII0=swi.Flux_FeIII0/((1-bsd.por)*bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
            
            swi.O20=100.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
            swi.NO30=40.0e-9;                                   % NO3 concentration at SWI (mol/cm^3)
            swi.NH40=10.0e-9;                                 	% NH4 concentration at SWI (mol/cm^3)
            swi.Fe20=0.0;                                       % Fe2 concentration at SWI (mol/cm^3)
            swi.SO40=2.8E-005;                                	% SO4 concentration at SWI (mol/cm^3)
            swi.H2S0=0.0;                                       % H2S concentration at SWI (mol/cm^3)
            swi.PO40=40.0e-9;                                   % PO4 concentration at SWI (mol/cm^3)
            swi.Mflux0=365*0.2e-10*1/(1-bsd.por)*1/bsd.w;       % actually CONCENTRATION of M at the sediment [mol/cm3] : from flux input  365*0.2e-10 (mol/(cm2*yr))
            swi.DIC0=2.4E-006;                                 	% DIC concentration at SWI (mol/cm^3)
            swi.ALK0=2.4E-006;                                 	% ALK concentration at SWI (mol/cm^3)
            swi.S0=35;                                         	% Salinity at SWI (not used at the moment)
            swi.plot_PO4_DIC_ALK=true;
        end
        
        function run_OMEN()
            % run OMEN-SED with default SWI conditions as in default_swi()
            clear
            tic;
            swi=benthic_test.default_swi()
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            toc;
            benthic_test.plot_column(res, false, res.swi, 'FULL_OMEN')
            
            % calculate depth integrated OM degradation rates
            Cox_rate.Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1, 1, res.bsd, res.swi, res);
            Cox_rate.Cox_Aerobic = res.zTOC.calcReac(0.0, res.zox, 1, 1, res.bsd, res.swi, res);
            if(res.swi.Nitrogen)
                Cox_rate.Cox_Denitr = res.zTOC.calcReac(res.zox, res.zno3, 1, 1, res.bsd, res.swi, res);
            end
            if(res.swi.Iron)
                Cox_rate.Cox_IronIII = res.zTOC.calcReac(res.zno3, res.zfeIII, 1, 1, res.bsd, res.swi, res);
            end
            Cox_rate.Cox_SO4red = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1, 1, res.bsd, res.swi, res)
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_10, C2_10] = res.zTOC.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC.calcOM(0.0, x, 1, 1, res.bsd, res.swi, res)
        end
        
        function run_PO4flux_SA_OMEN()
            % make PO4-SWI flux SA for changing boundary conditions in Corg and O2
            clear

            swi=benthic_test.default_swi()
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            for i=1:11   % 51
                swi.C01 = (0.1+(i-1)*0.2)*1e-2/12*2.5;     % 2.5 is rho_sed                      % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
                swi.C02 = swi.C01;
                Corg(i)=2*(0.1+(i-1)*0.2);
                for j=1:21  % 51              
                    swi.O20=(j-1)*3.0E-009;                                 % O2  concentration at SWI (mol/cm^3)                    
                    res=benthic_test.test_benthic(1,swi);
                    SWI_PO4(i,j)=res.flxswi_P;
                    O20(j)=(j-1)*3.0E-009;
%                  if((i==2 && j==2) || (i==3 && j==4))2.78E-005

%                      benthic_test.plot_column(res, false, swi, '00_expl_profiles')
%                  end
                end                
            end
            
            figure;
            hold on
            [C,h] = contourf(O20,Corg,SWI_PO4);
            clabel(C,h,'FontSize',16);
            colorbar();
            box on
            hold off
            xlabel('Ocean O_2')    
            ylabel ({'Corg'}); %;'(\mumol cm^{-2}yr^{-1})'})
            print('-depsc', '0_PO4_SWI-flux_SA_1811_NoAnoxicChange_FullP-cycle_higherPO4a_FasterDegrad_1000m');
%            print('-depsc', 'PO4_SWI-flux_SA_0907_4PO4_006500001_15Corg_150O2_3000m_ksfast');
        
        end
        
        function res = test_benthic( ncl, swi )

            if nargin < 1
                ncl = 1;
            end
            
            res.bsd = benthic_main(ncl);
            res.bsd.usescalarcode = ncl==1;
            
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            if ncl > 1  % set up O2 gradient for testing
                O20 = swi.O20;
                for i = 1:ncl
                    swi.O20(i) = 10*(i-1)/(ncl-1)*O20;
                end
            end
            
            res.swi = swi;
            
            % calculate
            res.zTOC = benthic_zTOC(res.bsd);
            res.zO2 = benthic_zO2(res.bsd, res.swi);
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            res.zFeIII = benthic_zFeIII(res.bsd, res.swi);
            res.zSO4 = benthic_zSO4(res.bsd, res.swi);
            res.zNH4 = benthic_zNH4(res.bsd, res.swi);
            res.zFe2 = benthic_zFe2(res.bsd, res.swi);
            res.zH2S = benthic_zH2S(res.bsd, res.swi);
            res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
            % with accelaration factor
%            res.zPO4_M = benthic_zPO4_M_Porg_CFA(res.bsd, res.swi);
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
            
            %            tic;
            % make zbio and k2 lower for anoxic conditions
         	if(swi.O20 < res.bsd.loc_BW_O2_anoxia)
%               res.bsd.zbio=0.0001;
%               res.zTOC.k2=0.00001; 
            end
            
            res = res.zTOC.calc(res.bsd,res.swi, res);
            % check O2 demand using O2 to C ratio and (convert POC concentr. to flux analog to fortran)
            % POC_flux*OC = POC_conc * w * 1/(1 - por) * OC
            O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por)
            
% %     change zbio and gammaH2S depending on BW oxygenation:
% %             if(res.swi.O20<=loc_BW_O2_anoxia)   % anoxic sediments?
% %                 bsd.zbio = 0.01;        % decrease bioturbation depth
% %                 bsd.gammaH2S = 0.95;    % fraction of H2S that is oxidised in anoxic sediments
% %             end
            
            if(res.swi.O20<=0.0)
                res.zox=0.0;
                res.flxzox = 0.0;
                res.conczox = 0.0;
                res.flxswiO2=0.0;
                res.zxf=0.0;
            else
                res = res.zO2.calc(res.bsd, res.swi, res);
            end
            if(swi.Nitrogen)
                res = res.zNO3.calc(res.bsd, res.swi, res);
            else
                res.zno3=res.zox;
                %                res.zso4=res.zox;   % for test-case with just TOC & O2
            end
            if(swi.Iron)
                res = res.zFeIII.calc(res.bsd, res.swi, res);
            else
                res.zfeIII = res.zno3;
            end
            res = res.zSO4.calc(res.bsd, res.swi, res);
            if(swi.Nitrogen)
                res = res.zNH4.calc(res.bsd, res.swi, res);
            end
            if(swi.Iron)
                                
                res = res.zFe2.calc(res.bsd, res.swi, res);
            end
            res = res.zH2S.calc(res.bsd, res.swi, res);
            res = res.zPO4_M.calc(res.bsd, res.swi, res);
            res = res.zDIC.calc(res.bsd, res.swi, res);
            res = res.zALK.calc(res.bsd, res.swi, res);
            %            toc;
            
            %%%%% WRITE OUTPUT:
            res.swi.C0nonbio = res.swi.C01nonbio+res.swi.C02nonbio;
            answ = res
            [Cinf, C1inf, C2inf] = res.zTOC.calcC( res.bsd.zinf, res.bsd, res.swi, res);
            [Cswi, C1swi, C2swi] = res.zTOC.calcC( 0, res.bsd, res.swi, res);
            fprintf('frac1 concentration at zinf %g \n',  C1inf);
            fprintf('frac2 concentration at zinf %g \n',  C2inf);
            fprintf('both concentration at zinf %g \n',  Cinf);
            fprintf(' \n');

            fprintf('frac1 concentration at swi nonbio %g \n',  res.swi.C01nonbio);
            fprintf('frac2 concentration at swi nonbio %g \n',  res.swi.C02nonbio);
            fprintf('both concentration at swi nonbio %g \n',  res.swi.C0nonbio);
            
            fprintf('sed preservation of POC %g \n',  Cinf/res.swi.C0nonbio);
            fprintf('DIC_SWI %.12g \n',  -(1-Cinf/res.swi.C0nonbio)*res.Fswi_TOC );
           
            fprintf(' \n');
            fprintf(' Calc with bioturbated [TOC] at SWI \n');
         	fprintf('frac1 concentration at swi bio %g \n',  C1swi);
            fprintf('frac2 concentration at swi bio %g \n',  C2swi);
            fprintf('both concentration at swi bio %g \n',  Cswi);
            fprintf('sed preservation of POC %g \n',  Cinf/Cswi);
            fprintf('wrong DIC_SWI %.12g \n',  -(1-Cinf/Cswi)*res.Fswi_TOC );
            fprintf(' \n');

            %             %%% WRITE EXACT FLUX
            %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
            %             fprintf('exact F_O2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
            
        end
        
        function plot_column(res, debug, swi, str_date)
            % plot single sediment column vs depth
            print_Fig1 = true;
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            % Nitrogen included?
            if(print_Fig1)
            if(swi.Nitrogen)
                %figure;
                fig1 = figure('Renderer', 'painters', 'Position', [10 10 600 900]);
                % TOC
                subplot(3,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                
                %%% TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                
                %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                %                ylim([-50 0.0])
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                %            title('Total TOC (wt%)')
                
                %%% O2
                if(res.zox>0.0)
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                else
                    O2(i) = 0.0;
                end
                subplot(3,2,3)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
%                ylim([-20 0.0])
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
                %%% NO3
                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,5)
                plot(NO3, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
%                ylim([-20 0.0])
                xlabel ('NO_3 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('NO3 (mol/cm^3)')
                
                %%% NH4
                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,2,4)
                plot(NH4, -zgrid, 'b')
                hold on
%                xlim([0 1e-7])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                xlabel ('NH_4 (mol/cm^3)')
                %                ylabel('Depth (cm)')
                %            title ('NH4 (mol/cm^3)')
                
                %%% SO4
                subplot(3,2,2)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                %                xlim([2.7e-5 swi.SO40])
                %                xlim([2.7e-5 swi.SO40])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                xlabel ('SO_4 (mol/cm^3)')
                %                ylabel('Depth (cm)')
                %            title ('SO4 (mol/cm^3)')
                
                %%% H2S
                subplot(3,2,6)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                %                xlim([0 4e-7])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                %               ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                % save Figure
                print(fig1, '-depsc2', ['0_ALL_PROFILES_' str_date '.eps']);
                
            else % no Nitrogen
                
                figure;
                % TOC
                subplot(2,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                %%% TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                
                %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                %            title('Total TOC (wt%)')
                
                %%% O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(2,2,3)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
                %%% SO4
                subplot(2,2,2)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
                %            xlim([2.7e-5 swi.SO40])
                xlabel ('SO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('SO4 (mol/cm^3)')
                
                %%% H2S
                subplot(2,2,4)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                
                print('-dpsc2', ['0_PROFILES_' str_date '.eps']);
            end
        end
            
            % Also plot the 2nd figure for PO4, DIC and ALK?
            if(swi.plot_PO4_DIC_ALK)
                % figure
            	fig2 = figure('Renderer', 'painters', 'Position', [10 10 600 900]);

                 %%% FeIII
                subplot(3,2,1)
                for i=1:length(zgrid)
                    [FeIII(i), flxFeIII(i)] = res.zFeIII.calcFeIII(zgrid(i), bsd, res.swi, res);
                end
                plot(FeIII, -zgrid, 'b')
                hold on
                %                xlim([2.7e-5 swi.FeIII0])
                %                xlim([2.7e-5 swi.FeIII0])
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                hold off
%                ylim_min = -round(res.zfeIII,-1)-10;
%                ylim([ylim_min 0.0])
                xlabel ('FeIII (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('SO4 (mol/cm^3)')

                %%% Fe2
                subplot(3,2,2)
                for i=1:length(zgrid)
                    [Fe2(i), flxFe2(i)] = res.zFe2.calcFe2(zgrid(i), bsd, res.swi, res);
                end
                plot(Fe2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('Fe2 (mol/cm^3)')
%                ylabel('Depth (cm)')
                %            title ('Fe2 (mol/cm^3)')
                
                
                %%% PO4
                subplot(3,2,3)
                for i=1:length(zgrid)
                    [PO4(i), flxPO4(i), M(i), flxM(i), e_M(i), f_M(i), p_M(i), q_M(i), g_M(i), dedz_M(i), dfdz_M(i), dpdz_M(i), dqdz_M(i), dgdz_M(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                end
                plot(PO4, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                %          axis([0 1.5*10^(-9) -100 0])
                xlabel ('PO_4 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('PO_4 (mol/cm^3)')
                
                %%% Fe-bound P (M)
                subplot(3,2,4)
                %for i=1:length(zgrid)
                %    [PO4(i), flxPO4(i), M(i), flxM(i)] = res.zPO4_M.calcPO4_M(zgrid(i), bsd, res.swi, res);
                %end
                plot(M, -zgrid, 'b')
                hold on
                %            plot([0,max(M)], [-bsd.zbio,-bsd.zbio], 'k--')
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('Fe-bound P (mol/cm^3)')
                %            ylabel('Depth (cm)')
                %            title ('Fe-bound P (mol/cm^3)')
                
                %%% DIC
                subplot(3,2,5)
                for i=1:length(zgrid)
                    [DIC(i), flxDIC(i)] = res.zDIC.calcDIC(zgrid(i), bsd, res.swi, res);
                end
                plot(DIC, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('DIC (mol/cm^3)')
                ylabel('Depth (cm)')
                
                %%% ALK
                subplot(3,2,6)
                for i=1:length(zgrid)
                    [ALK(i), flxALK(i)] = res.zALK.calcALK(zgrid(i), bsd, res.swi, res);
                end
                plot(ALK, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('ALK (mol/cm^3)')
                ylabel('Depth (cm)')
                
                
                print(fig2, '-depsc2', ['0_PO4_PROFILES_' str_date '.eps']);
            end
            
            
            %%% some plots for debugging!
            if debug
                figure
                
                subplot(3,2,1)
                hold on
                plot(e_M, -zgrid, 'b')
                title ('Fe-P - ODE solution')
                
                subplot(3,2,2)
                hold on
                plot(f_M, -zgrid, 'b')
                
                subplot(3,2,3)
                hold on
                plot(p_M, -zgrid, 'b')
                
                subplot(3,2,4)
                hold on
                plot(q_M, -zgrid, 'b')
                
                subplot(3,2,5)
                hold on
                plot(g_M, -zgrid, 'b')
                
                figure
                hold on
                subplot(3,2,1)
                plot(dedz_M, -zgrid, 'b')
                title ('Fe-P - ODE derivations')
                subplot(3,2,2)
                plot(dfdz_M, -zgrid, 'b')
                subplot(3,2,3)
                plot(dpdz_M, -zgrid, 'b')
                subplot(3,2,4)
                plot(dqdz_M, -zgrid, 'b')
                subplot(3,2,5)
                plot(dgdz_M, -zgrid, 'b')
                
                
                %%%%%%%%%%%%%%%%%%%%%
                
                %         H2S
                
                %%%%%%%%%%%%%%%%%%%%%
                figure
                subplot(3,3,2)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S_debug(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                xlabel ('H_2S (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('H2S (mol/cm^3)')
                
                subplot(3,3,4)
                hold on
                plot(e_H2S, -zgrid, 'b')
                title ('H_2S - ODE solution')
                
                subplot(3,3,5)
                hold on
                plot(f_H2S, -zgrid, 'b')
                
                subplot(3,3,6)
                hold on
                plot(g_H2S, -zgrid, 'b')
                
                subplot(3,3,7)
                hold on
                plot(dedz_H2S, -zgrid, 'b')
                title ('H_2S - ODE derivations')
                
                subplot(3,3,8)
                hold on
                plot(dfdz_H2S, -zgrid, 'b')
                
                subplot(3,3,9)
                hold on
                plot(dgdz_H2S, -zgrid, 'b')
                
                
                
                % CONCENTRATION + Vertical Trransport
                figure;
                % TOC
                subplot(3,4,1)
                for i=1:length(zgrid)
                    [C(i), C1(i), C2(i)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i), C2flx(i)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                % TOC wt %
                plot(100*C1*12/bsd.rho_sed, -zgrid, 'b')
                hold on
                plot(100*C2*12/bsd.rho_sed, -zgrid, 'g')
                plot(100*C*12/bsd.rho_sed, -zgrid, 'k')
                plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                title('Total TOC (wt%)')
                % TOC vertical transport flux
                subplot(3,4,2);
                plot(C1flx, -zgrid, 'b')
                hold on
                plot(C2flx, -zgrid, 'g')
                plot(Cflx, -zgrid, 'k')
                xlabel ('TOC trspt (mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title('TOC vert transport')
                
                
                % O2
                for i=1:length(zgrid)
                    [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,3)
                plot(O2, -zgrid, 'b')
                hold on
                plot([0,res.swi.O20], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('O2 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('O2 (mol/cm^3)')
                subplot(3,4,4);
                plot(flxO2, -zgrid, 'b', flxO2D,-zgrid,'b--',flxO2adv,-zgrid,'c--');%,flxO2D+flxO2adv,-zgrid,'r--');
                legend('tot','diff','adv','diff+adv');
                legend boxoff;
                xlabel ('O2 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('O2 vert transport')
                
                
                %             % NO3
                %
                %             for i=1:length(zgrid)
                %                 [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                %             end
                %             subplot(3,4,5)
                %             plot(NO3, -zgrid, 'b')
                %             hold on
                %             plot([0,res.swi.NO30], [-bsd.zbio,-bsd.zbio], 'k--')
                %             xlabel ('NO3 (mol/cm^3)')
                %             ylabel('Depth (cm)')
                %             title ('NO3 (mol/cm^3)')
                %             subplot(3,4,6)
                %             plot(flxNO3, -zgrid, 'b')
                %             xlabel ('NO3 trsp(mol cm^{-2}yr^{-1})')
                %             ylabel('Depth (cm)')
                %             title ('NO3 vert transport');
                %
                %
                %
                %             for i=1:length(zgrid)
                %                 [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                %             end
                %             subplot(3,4,7)
                %             plot(NH4, -zgrid, 'b')
                %             hold on
                %             plot([0,res.swi.NH40], [-bsd.zbio,-bsd.zbio], 'k--')
                %             hold off
                %             xlabel ('NH4 (mol/cm^3)')
                %             ylabel('Depth (cm)')
                %             title ('NH4 (mol/cm^3)')
                %             subplot(3,4,8)
                %             plot(flxNH4, -zgrid, 'b');
                %             xlabel ('NH4 trsp(mol cm^{-2}yr^{-1})')
                %             ylabel('Depth (cm)')
                %             title ('NH4 vert transport')
                
                subplot(3,4,9)
                for i=1:length(zgrid)
                    [SO4(i), flxSO4(i)] = res.zSO4.calcSO4(zgrid(i), bsd, res.swi, res);
                end
                plot(SO4, -zgrid, 'b')
                hold on
                plot([0,res.swi.SO40], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                %xlim([0 SO40])
                xlabel ('SO4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('SO4 (mol/cm^3)')
                subplot(3,4,10)
                plot(flxSO4, -zgrid, 'b');
                xlabel ('SO4 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('SO4 vert transport')
                
                subplot(3,4,11)
                for i=1:length(zgrid)
                    [H2S(i), flxH2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                end
                plot(H2S, -zgrid, 'b')
                hold on
                plot([0,res.swi.H2S0], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('H2S (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('H2S (mol/cm^3)')
                subplot(3,4,12)
                plot(flxH2S, -zgrid, 'b');
                xlabel ('H2S trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('H2S vert transport')
                
            end
            
        end
        
    end
    
end

