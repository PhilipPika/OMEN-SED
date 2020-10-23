%% OMEN-SED 1.1 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al. (2017) GMD
% including the RCM approximation of Pika et al. (2020) GMD

% benthic_test.m
% functions to run OMEN-SED and plot the results

% Command to run the model: 
% benthic_test.run_OMEN: For original 2G model
% benthic_test.run_OMEN_RCM: For RCM approximation
% 1) default sediment-water interface boundary conditions are set as prescribed in default_swi()
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
            
            swi.Test_Dale = true;
            swi.plot_fig = false;                                % plot the sediment profiles
            
            swi.Nitrogen=true;                                  % calculate N (true/false)
            swi.Iron=true;                                      % calculate Fe (true/false)
            % for 2G-model
            swi.C01_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 solid phase)
            swi.C02_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 solid phase)
            swi.Fnonbio1 =  0.5*6.2e-3/100^2*365;    % swi.C01_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.Fnonbio2 =  0.5*6.2e-3/100^2*365;    % swi.C02_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.C01 = swi.C01_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.C02 = swi.C02_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m           
            
            % for nG-model
            swi.nG = 100;
            swi.p_a =100.1;  % 3e-4;
            swi.p_nu = 0.125;
            swi.C0_nonbio = 6.0 * 1e-2/12*bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)


%            swi.FeIII0=2.8E-005; %3.0E-006;                   	% FeIII concentration at SWI (mol/cm^3) --> TODO: needs to be a flux!
            Fe_influx_Dale = 0.1667*1110.0E-006;                     % FeIII influx from Dale but just 16.67% is available for DIR (See his Tab. 2, footnote d)
            swi.Flux_FeIII0 =  (Fe_influx_Dale)*365/100^2;           % Dale 1110 mol/(m^2 day)   -->  mol/(cm^2 yr):     *365/100^2
%            swi.FeIII0= swi.Flux_FeIII0/((1-bsd.por)*bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
            
            swi.O20=200.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
            swi.NO30=35.0e-9;                                   % NO3 concentration at SWI (mol/cm^3)
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
            swi=benthic_test.default_swi();
            swi.TwoG_OM_model = true;
            %            % set date-time
            %            str_date = datestr(now,'ddmmyy_HH_MM_SS');
            res=benthic_test.test_benthic(1,swi);
            toc;
            benthic_test.plot_column(res, false, res.swi, 'FULL_OMEN')
            
            % calculate depth integrated OM degradation rates (wrong use
            % difference of fluxes (see test_benthic() below)
            Cox_rate.Cox_total = res.zTOC.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            Cox_rate.Cox_Aerobic = res.zTOC.calcReac(0.0, res.zox, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            if(res.swi.Nitrogen)
                Cox_rate.Cox_Denitr = res.zTOC.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            end
            if(res.swi.Iron)
                Cox_rate.Cox_IronIII = res.zTOC.calcReac(res.zno3, res.zfeIII, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res);
            end
            Cox_rate.Cox_sulfate = res.zTOC.calcReac(res.zno3, res.bsd.zinf, 1*(1-res.bsd.por), 1*(1-res.bsd.por), res.bsd, res.swi, res)
            
            % Calculate out put to compare with Dale:
%            Cox_total_Dale_units = Cox_rate.Cox_total*1000 *100^2/365      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a

            % calculate mean OM concentration in upper x cm
            [C_10, C1_10, C2_10] = res.zTOC.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC.calcOM(0.0, x, 1, 1, res.bsd, res.swi, res);
        end
        
        function [Cox_rate, Flux_Fe2_Dale_units, res, Output]  = run_OMEN_RCM()
            % run OMEN-SED with default SWI conditions as in default_swi()
            clear
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            
            if(swi.Test_Dale)
                swi.nG = 14;
                swi.p_a = 0.1; 
                swi.p_nu = 0.125;
              	swi.POC_flux = [0.5 1 2 4 6 8 10 12 14 16];
              	swi.POCi = 5;                
            end

            res=benthic_test.test_benthic(1,swi);
            % set date-time or string going into plot function
            str_date = [num2str(res.swi.nG) 'G_a=' num2str(res.swi.p_a) '_nu=' num2str(res.swi.p_nu) '_O20_' num2str(res.swi.O20)];
            if(swi.plot_fig)
                benthic_test.plot_column(res, false, swi, str_date)
                benthic_test.plot_TOC(res, false, swi, str_date)
            end
            % calculate depth integrated OM degradation rates (this is
            % wrong - use the fluxes of TOC at the boundaries - see test_benthic)
            swi.C0i = res.swi.C0i;
            Cox_rate(1,1) = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res);
            Cox_rate(1,2) = res.zTOC_RCM.calcReac(0.0, res.zox, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            if(swi.Nitrogen)
                Cox_rate(1,3) = res.zTOC_RCM.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            end
            if(res.swi.Iron)
                Cox_rate(1,4) = res.zTOC_RCM.calcReac(res.zno3, res.zfeIII,1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;
            end
            Cox_rate(1,5) = res.zTOC_RCM.calcReac(res.zfeIII,res.zso4, 1*(1-res.bsd.por), res.bsd, swi, res)/Cox_rate(1,1)*100;

            if(swi.Test_Dale)
                % compare results with Dale
                [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);

                Cox_total = (F_TOC_swi-F_TOC_inf)*1000 *100^2/365;      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
                Flux_Fe2_Dale_units = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                fprintf('\n');
                fprintf('Total Cox as flux difference (mmol m-2 d-1) %e \n',  Cox_total);
                fprintf('SWI-flux FeII (umol m-2 d-1) %e \n',  Flux_Fe2_Dale_units);

            end
            
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_11] = res.zTOC_RCM.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, x, 1, res.bsd, swi, res);
            
            Output(1,1) = res.flxswiH2S;
            Output(1,2) = res.flxswiSO4;
            Output(1,3) = res.flxswiFe2;
            Output(1,4)  = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
            Output(1,5) = res.flxswiFeIII;
            Output(1,6) = res.flxswiNH4;
            Output(1,7) = res.flxswiNO3;
            Output(1,8) = res.flxswiO2;
            Output(1,9) = res.zox;
            Output(1,10) = res.zno3;
            Output(1,11) = res.zfeIII;
            Output(1,12) = res.zso4;
        end
        
        
        function [Cox_rate_out, Penetration_out, Flux_Fe2_Dale_units] = run_OMEN_RCM_testDale_POC()
            % run OMEN-SED with default SWI conditions as in default_swi()
            clear
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            
            if(swi.Test_Dale)
                swi.nG = 14;
                swi.p_a = 0.1; 
                swi.p_nu = 0.125;
                
            end
            
            swi.POC_flux = [0.5 1 2 4 6 8 10 12 14 16];

            for i=1:length(swi.POC_flux)
                swi.POCi = i;
                res=benthic_test.test_benthic(1,swi);
                % set date-time or string going into plot function

                swi.C0i = res.swi.C0i;
                Cox_rate.Cox_aerobic = res.zTOC_RCM.calcReac(0.0, res.zox, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                if(swi.Nitrogen)
                    Cox_rate.Cox_denitr = res.zTOC_RCM.calcReac(res.zox, res.zno3, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                end
                if(res.swi.Iron)
                    Cox_rate.Cox_FeIII = res.zTOC_RCM.calcReac(res.zno3, res.zfeIII,1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                end
                Cox_rate.Cox_sulfate = res.zTOC_RCM.calcReac(res.zfeIII,res.zso4, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                Cox_rate.Cox_total = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365;
                Cox_rate.Total_rate_sum = Cox_rate.Cox_aerobic + Cox_rate.Cox_denitr + Cox_rate.Cox_FeIII + Cox_rate.Cox_sulfate;
                % calculate degradation pathways
                Cox_rate_out(1,2*i-1) = Cox_rate.Cox_aerobic;
                Cox_rate_out(2,2*i-1) = Cox_rate.Cox_denitr;
                Cox_rate_out(3,2*i-1) = Cox_rate.Cox_FeIII;
                Cox_rate_out(4,2*i-1) = Cox_rate.Cox_sulfate;
                Cox_rate_out(5,2*i-1) = Cox_rate.Cox_total;
                % FRaction of total:
                Cox_rate_out(1,2*i) = Cox_rate.Cox_aerobic/Cox_rate.Cox_total*100;
                Cox_rate_out(2,2*i) = Cox_rate.Cox_denitr/Cox_rate.Cox_total*100;
                Cox_rate_out(3,2*i) = Cox_rate.Cox_FeIII/Cox_rate.Cox_total*100;
                Cox_rate_out(4,2*i) = Cox_rate.Cox_sulfate/Cox_rate.Cox_total*100;
                Cox_rate_out(5,2*i) = Cox_rate.Cox_total/Cox_rate.Cox_total*100;


                Penetration_out(1,2*i-1) = res.zox;
                Penetration_out(2,2*i-1) = res.zno3;
                Penetration_out(3,2*i-1) = res.zfeIII;
                Penetration_out(4,2*i-1) = res.zso4;
                Penetration_out(1,2*i) = nan;
                Penetration_out(2,2*i) = nan;
                Penetration_out(3,2*i) = nan;
                Penetration_out(4,2*i) = nan;
                % compare results with Dale
                [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);

                Cox_total_diff = (F_TOC_swi-F_TOC_inf)*1000 *100^2/365;      % in units of mmol m-2 d-1 as in Dale ea. 2015, Fig. 2a
                Flux_Fe2_Dale_units(2*i-1) = res.flxswiFe2 *10^6 *100^2 /365;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
                Flux_Fe2_Dale_units(2*i) = nan;                   % in units of umol m-2 d-1 as in Dale ea. 2015, Fig. 3
            end
        end
        
        
        function res = test_benthic( ncl, swi )
%            loc_BW_O2_anoxia = 5.0e-9;       	% set to 5.0 nanomol/cm^3
            if nargin < 1
                ncl = 1;
            end
            
            if(swi.Test_Dale)
                res.bsd = benthic_main(ncl, 350);   % specify 350m water-depth
            else
                res.bsd = benthic_main(ncl);                
            end
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
            
            % initialize then calculate
            if(swi.TwoG_OM_model)
                res.zTOC = benthic_zTOC(res.bsd);
            else
                res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
            end
   
            tic;
            if(swi.TwoG_OM_model)
                res = res.zTOC.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por);
            else
                % Adding into on RCM for MultiG approach
                if(res.swi.Test_Dale)
                    [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM_Dale_2015(res.bsd, res.swi);
                else
                    [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                end
                res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(sum(res.swi.Fnonbioi))*res.bsd.OC/((1-res.bsd.por)./res.bsd.por);
            end
            Cox_total = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1*(1-res.bsd.por), res.bsd, swi, res)*1000 *100^2/365; % calculate in mmol/(m2 d) for Seb's parameteriztion
            res.zO2 = benthic_zO2(res.bsd, res.swi);
            res.zNO3 = benthic_zNO3(res.bsd, res.swi);
            res.zFeIII = benthic_zFeIII(res.bsd, res.swi);                      
            res.zSO4 = benthic_zSO4(res.bsd, res.swi);
            res.zNH4 = benthic_zNH4(res.bsd, res.swi);
            res.zFe2 = benthic_zFe2(res.bsd, res.swi, Cox_total);
            res.zH2S = benthic_zH2S(res.bsd, res.swi);
            res.zPO4_M = benthic_zPO4_M(res.bsd, res.swi);
            res.zDIC = benthic_zDIC(res.bsd, res.swi);
            res.zALK = benthic_zALK(res.bsd, res.swi);
         
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
                % First solve iteratively for correct SWI [Fe3+]
                fun=@(Fe3_conc_in) -res.zFeIII.calc_input(res.bsd, res.swi, res, Fe3_conc_in) - res.swi.Flux_FeIII0;   % NOTE: negative calc_input to transform the result to positive because the given flux at the SWI is positive here!
                
                % For starting point calculate SWI-concentration of Fe3+ from input flux
                if(res.zno3 > res.bsd.zbio )
                    % SWI-concentration of Fe3+ is not affected by Fe-reduction
                    % -> use simple conversion - this should be close
                    % enough to correct SWI-concentration
                    res.swi.FeIII0= res.swi.Flux_FeIII0/((1-res.bsd.por)*res.bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
                    test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                 	% try zero flux at zinf and see if we have any FeIII left, also
                    % calculate [FeIII] at zinf for advective loss
                    [flxzfeIII, conczinf, flxswi,rtmp] = res.zFeIII.calcbc(res.bsd.zinf, res.bsd, res.swi, res, 2); 
                    if(conczinf >= 0)    % zfeIII >= zinf --> use first guess for fe3-flux and assume the remainder is buried                        
                        fun0 = 0.0;     % use the first guess for res.swi.FeIII0
                        % FOR MASS BALANCE WITH cGENIE MAKE SURE THE
                        % REMAINING INFLUX OF Fe3+ IS BURIED IN SEDIMENTS!!!
                    else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
                        fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration     
                    end
                else
                    % SWI-concentration of Fe3+ is affected by Fe-reduction
                    % -> take biodiffusion and Fe-reduction into account
                    res.swi.FeIII0 = res.zFeIII.calcConcentration(res.bsd, res.swi, res, res.swi.Flux_FeIII0);
                    if(res.swi.FeIII0 < 0)  
                        % Case zno3 < zFe3 < zbio: calcConcentration() fails
                        % use simple conversion and iterate from there
                        res.swi.FeIII0= res.swi.Flux_FeIII0/((1-res.bsd.por)*res.bsd.w);     % calculate concentration [mol/cm^3] from flux [mol/(cm2 yr)] according non-bioturbated flux!!!
                         test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                        fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration                    
                    else  
                 	% try zero flux at zinf and see if we have any FeIII left, also
                    % calculate [FeIII] at zinf for advective loss
                    [flxzfeIII, conczinf, flxswi,rtmp] = res.zFeIII.calcbc(res.bsd.zinf, res.bsd, res.swi, res, 2); 
                    if(conczinf >= 0)    % zfeIII >= zinf --> use first guess for fe3-flux and assume the remainder is buried                        
                        fun0 = 0.0;     % use the first guess for res.swi.FeIII0
                    else    % zfeIII < zinf --> iteratively solve for correct Fe3-concentration
%                         test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                        fun0 = fun(res.swi.FeIII0);     % initial test with calculated SWI-concentration     
                    end
                    end
                end
                
                if(abs(fun0) > res.bsd.tol_Fe3)  % initial guess above not good enough
                    
                    if(fun0 > 0)    % Fe3+ concentration too high
%                         test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 1e-10);
%                         test1_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);
                        res.swi.FeIII0=fzero(fun,[1e-10, res.swi.FeIII0],res.bsd.fzerooptions_Fe3);    % NOTE: Find better way for lower boundary
                    else    % Fe3+ concentration too low
%                          test0_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, res.swi.FeIII0);                        
%                          test1_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 0.5*res.swi.FeIII0);
%                          test2_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 10*res.swi.FeIII0);
%                          test3_flx = res.zFeIII.calc_input(res.bsd, res.swi, res, 1e+6);
                        res.swi.FeIII0=fzero(fun,[res.swi.FeIII0, 1e+6],res.bsd.fzerooptions_Fe3);    % NOTE: Find better way for upper boundary - not sure if this works anyway (at least not for gammaFe2 >= 0.99)
                    end                    
                end                    
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
            toc;
            
            %%%%% WRITE OUTPUT:
            answ = res
            if(swi.TwoG_OM_model)
                % use fluxes to calculate burial efficiency and total Cox:
                [F_TOC_swi, F_TOC1_swi] = res.zTOC.calcCflx(0, res.bsd, res.swi, res);
                [F_TOC_inf, F_TOC1_inf] = res.zTOC.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);

                [F_TOC_zox, F_TOC1_zox] = res.zTOC.calcCflx(res.zox, res.bsd, res.swi, res);
                [F_TOC_zno3, F_TOC1_zno3] = res.zTOC.calcCflx(res.zno3, res.bsd, res.swi, res);
                [F_TOC_zfe3, F_TOC1_zfe3] = res.zTOC.calcCflx(res.zfeIII, res.bsd, res.swi, res);
                [F_TOC_zso4, F_TOC1_zso4] = res.zTOC.calcCflx(res.zso4, res.bsd, res.swi, res);

                fprintf('Fluxes at zinf %g \n',  F_TOC_inf);
                fprintf('Fluxes at swi %g \n',  F_TOC_swi);
                
                fprintf('Burial efficiency of POC (in %%) %g \n',  F_TOC_inf/F_TOC_swi *100);
                fprintf('\n');
                fprintf('Total Cox (mol cm-2 yr-1) %e \n',  F_TOC_swi-F_TOC_inf);
                fprintf('Aerobic Cox %e \n',  F_TOC_swi-F_TOC_zox);
                fprintf('Denitrification Cox %e \n',  F_TOC_zox - F_TOC_zno3 );
                fprintf('Fe-reduction Cox %e \n',  F_TOC_zno3 - F_TOC_zfe3);
                fprintf('Sulfate reduction Cox %e \n',  F_TOC_zfe3-F_TOC_zso4);

            else
                % use fluxes to calculate burial efficiency and total Cox:
                [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(res.bsd.zinf, res.bsd, res.swi, res);
                
                [F_TOC_zox, F_TOC1_zox] = res.zTOC_RCM.calcCflx(res.zox, res.bsd, res.swi, res);
                [F_TOC_zno3, F_TOC1_zno3] = res.zTOC_RCM.calcCflx(res.zno3, res.bsd, res.swi, res);
                [F_TOC_zfe3, F_TOC1_zfe3] = res.zTOC_RCM.calcCflx(res.zfeIII, res.bsd, res.swi, res);
                [F_TOC_zso4, F_TOC1_zso4] = res.zTOC_RCM.calcCflx(res.zso4, res.bsd, res.swi, res);

                fprintf('Fluxes at zinf %g \n',  F_TOC_inf);
                fprintf('Fluxes at swi %g \n',  F_TOC_swi);                                
                fprintf('Burial efficiency of POC (in %%) %g \n',  F_TOC_inf/F_TOC_swi *100);
                fprintf('\n');
                fprintf('Total Cox (mol cm-2 yr-1) %e \n',  F_TOC_swi-F_TOC_inf);
                fprintf('Aerobic Cox %e \n',  F_TOC_swi-F_TOC_zox);
                fprintf('Denitrification Cox %e \n',  F_TOC_zox - F_TOC_zno3 );
                fprintf('Fe-reduction Cox %e \n',  F_TOC_zno3 - F_TOC_zfe3);
                fprintf('Sulfate reduction Cox %e \n',  F_TOC_zfe3-F_TOC_zso4);
                
                %             %%% WRITE EXACT FLUX
                %             FO2_exact=res.zO2.calcFO2_exact(res.zox,res.bsd, res.swi, res);
                %             fprintf('exact flxswiO2 flux (mol cm^{-2} yr^{-1}) %g \n',  FO2_exact);
                
            end
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
                if(swi.TwoG_OM_model)
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
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    %            plot([0,(res.swi.C01+res.swi.C02)*12/bsd.rho_sed ], [-bsd.zbio,-bsd.zbio], 'k--')
                    hold off
                    %                ylim([-50 0.0])
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    %            title('Total TOC (wt%)')
                else
                    subplot(3,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    %                 color = parula(swi.nG);
                    %                 for G = 1:swi.nG
                    %                     plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                    %                     hold on
                    %                 end
                    % Plot sum (TOC)
                    plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    %            title('Total TOC (wt%)')
                    
                end
                
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
                if(swi.TwoG_OM_model)
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
                else
                    subplot(2,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    color = parula(swi.nG);
                    for G = 1:swi.nG
                        plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        hold on
                    end
                    % Plot sum (TOC)
                    plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    %            title('Total TOC (wt%)')
                    
                end
                
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
                    [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S(zgrid(i), bsd, res.swi, res);
                    % no func. called calcH2S_debug [H2S(i), flxH2S(i), e_H2S(i), dedz_H2S(i), f_H2S(i), dfdz_H2S(i), g_H2S(i), dgdz_H2S(i)] = res.zH2S.calcH2S_debug(zgrid(i), bsd, res.swi, res);
                    
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
                if(swi.TwoG_OM_model)
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
                else
                    subplot(3,4,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    color = parula(swi.nG);
                    for G = 1:swi.nG
                        plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        hold on
                    end
                    % Plot sum (TOC)
                    plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                    hold on
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    % TOC vertical transport flux
                    subplot(3,4,2);
                    plot(C1flx, -zgrid, 'b')
                    hold on
                    plot(Cflx, -zgrid, 'k')
                    xlabel ('TOC trspt (mol cm^{-2}yr^{-1})')
                    ylabel('Depth (cm)')
                    title('TOC vert transport')
                    
                end
                
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
                
                
                for i=1:length(zgrid)
                    [NO3(i), flxNO3(i)] = res.zNO3.calcNO3(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,5)
                plot(NO3, -zgrid, 'b')
                hold on
                plot([0,res.swi.NO30], [-bsd.zbio,-bsd.zbio], 'k--')
                xlabel ('NO3 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('NO3 (mol/cm^3)')
                subplot(3,4,6)
                plot(flxNO3, -zgrid, 'b')
                xlabel ('NO3 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('NO3 vert transport');
                
                
                
                for i=1:length(zgrid)
                    [NH4(i), flxNH4(i)] = res.zNH4.calcNH4(zgrid(i), bsd, res.swi, res);
                end
                subplot(3,4,7)
                plot(NH4, -zgrid, 'b')
                hold on
                plot([0,res.swi.NH40], [-bsd.zbio,-bsd.zbio], 'k--')
                hold off
                xlabel ('NH4 (mol/cm^3)')
                ylabel('Depth (cm)')
                title ('NH4 (mol/cm^3)')
                subplot(3,4,8)
                plot(flxNH4, -zgrid, 'b');
                xlabel ('NH4 trsp(mol cm^{-2}yr^{-1})')
                ylabel('Depth (cm)')
                title ('NH4 vert transport')
                % Till here
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
                
                
                if(~swi.TwoG_OM_model)
                    figure;
                    %%% TOC
                    subplot(3,2,1)
                    for i=1:length(zgrid)
                        [C(i), C1(i,:)] = res.zTOC.calcC( zgrid(i), bsd, res.swi, res);
                        [Cflx(i), C1flx(i,:)] = res.zTOC.calcCflx( zgrid(i), bsd, res.swi, res);
                    end
                    color = parula(swi.nG);
                    for G = 1:swi.nG
                        plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                        hold on
                    end
                    
                    t=xlim;         % to draw penetration depths the correct lengths
                    plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                    plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                    plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
                    plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
                    
                    hold off
                    xlabel ('TOC (wt%)')
                    ylabel('Depth (cm)')
                    title('Total TOC (wt%)')
                end
            end
            
        end
        
        function plot_TOC(res, ~, swi, str_date)
            figure
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            %%% TOC
            for i=1:length(zgrid)
                [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            subplot(121)
            color = parula(swi.nG);
            for G = 1:swi.nG
                plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                hold on
            end
            % Plot sum (TOC)
            subplot(122)
            plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
            hold on
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
            plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
            plot([0,t(1,2)], [-res.zno3,-res.zno3], 'g--')
            plot([0,t(1,2)], [-res.zfeIII,-res.zfeIII], 'm--')
            plot([0,t(1,2)], [-res.zso4,-res.zso4], 'r--')
            
            hold off
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            %            title('Total TOC (wt%)')
        end
        
        function plot_TOC_O2_column(res, debug, swi, str_date)
            % plot single sediment column vs depth
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:0.1:bsd.zinf;
            
            % CONCENTRATIONS WITHOUT PO4
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            % Nitrogen included?
            if(swi.Nitrogen)
                figure;
                
                %%% TOC
                subplot(1,2,1)
                for i=1:length(zgrid)
                    [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                    [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
                end
                color = parula(swi.nG);
                for G = 1:swi.nG
                    plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                    hold on
                end
                % Plot sum (TOC)
                plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                
                hold off
                xlabel ('TOC (wt%)')
                ylabel('Depth (cm)')
                
                
                text( 3/5*(100*C(1)*12/bsd.rho_sed) ,-60,{['a:' num2str(swi.p_a)],...
                    ['SFD:' num2str(bsd.wdepth)]})
                
                %%% O2
                if(res.zox>0.0)
                    for i=1:length(zgrid)
                        [O2(i), flxO2(i), flxO2D(i), flxO2adv(i)] = res.zO2.calcO2(zgrid(i), bsd, res.swi, res);
                    end
                else
                    O2(i) = 0.0;
                end
                subplot(1,2,2)
                plot(O2, -zgrid, 'b')
                hold on
                t=xlim;         % to draw penetration depths the correct lengths
                plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
                plot([0,t(1,2)], [-res.zox,-res.zox], 'b--')
                %                 ylim([-20 0.0])
                xlabel ('O_2 (mol/cm^3)')
                ylabel('Depth (cm)')
                %            title ('O2 (mol/cm^3)')
                
            end
        end       
                
        function [k, C0i, Fnonbioi] = RCM(bsd, swi)
            %% Implemented by Philip Pika
            %% Described in Pika et al. (2020) GMD
            
            % For comparison with original 2G results
            if swi.nG == 2          % nothing needs to be done here - all set in default_swi()
%                 F(1) = 0.5;
%                 F(2) = 0.5;
%                 C0i = F.*swi.C0_nonbio * 1e-2/12*bsd.rho_sed;       % TOC@SWI (wt%) -> (mol/cm^3 bulk phase), 2.5 sed.density (g/cm3) 0.1
%                 Fnonbioi = swi.C0_nonbio*(1-bsd.por)*bsd.w;        % [mol/(cm2 yr)] according non-bioturbated flux
%                 k = [0.01 0.0001];
%                 swi.p_a = NaN;
%                 swi.p_nu = NaN;
            else
                emin = -15; % as in Dale ea. '15: -10;
                emax = -log10(swi.p_a)+2;  % as in Dale ea. '15: 2    % upper k limit for k-bins of multi-G approximation
                if emin >= emax;error('emin >= emax, this cannot be!');end
%                 if emax >= log10(200);emax=log10(200);end
                
                k(1)= 10^(emin);
                kk(1)=10^(emin);
                F(1) = gammainc(swi.p_a*10^emin,swi.p_nu,'lower');
                kk(swi.nG)=10^(emax);
                k(swi.nG)=10^(emax);
                F(swi.nG) = gammainc(swi.p_a*10^emax,swi.p_nu,'upper');
                
                % Define the b.c. for all the intermediate fractions
                
                G=2:swi.nG-1;
                ne=emin+(1:swi.nG-2).*(emax-emin)./(swi.nG-1);
                kk(2:swi.nG-1)=10.^ne;
                G_inc_0 = gammainc(swi.p_a*kk(1:swi.nG-2),swi.p_nu,'upper'); % G-1 = 1:end-2
                G_inc_1 = gammainc(swi.p_a*kk(2:swi.nG-1),swi.p_nu,'upper'); % G = 2:end-1
                F(2:swi.nG-1) = (G_inc_0 - G_inc_1);
                k(G)=kk(1:swi.nG-2)+(kk(2:swi.nG-1)-kk(1:swi.nG-2))/2;
                F(F<=eps)=eps;
                if abs(sum(F)-1) > 0.0001;warning('F~=1!!');end
                Fnonbioi = F.* ( swi.C0_nonbio*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
                C0i = F.*swi.C0_nonbio;
            end
        end
        
      	function [k, C0i, Fnonbioi] = RCM_Dale_2015(bsd, swi)
            %% specify k and F speficically as in Dale ea. 2015 and use input flux to calculate non-bioturbated SWI-concentration of POC
            

                k(1)= 1e-10;
                for i=2:13
                    j = -12+i;
                    k(i)= 3.16*10^j;
                end
                k(14)= 100;
%                k(1:10) = 0.08;    % in case we want to degrade all
                
                F(1) = 0.021746;
                F(2) = 0.00725275;
                F(3) = 0.0096717;
                F(4) = 0.0128974;
                F(5) = 0.017199;
                F(6) = 0.0229352;
                F(7) = 0.0305846;
                F(8) = 0.0407852;
                F(9) = 0.0543879;
                F(10) = 0.0725265;
                F(11) = 0.0967046;
                F(12) = 0.12881;
                F(13) = 0.169822;
                F(14) = 0.314677;
                
                if abs(sum(F)-1) > 0.0001;warning('F~=1!!');end
                
                % use the input flux of 6.2 mmol m-2 d-1 from Dale
                % BC Dale ea 2015: 
%                POC_flux = [0.5 1 2 4 6 8 10 12 14 16];

                FPOC_Dale = swi.POC_flux(swi.POCi)*10^(-3)/100^2*365;   %   in mol / (cm^2*yr)
%                FPOC_Dale = 6.0*10^(-3)/100^2*365;   %   in mol / (cm^2*yr)
                Fnonbioi = F.* FPOC_Dale; % NonBioturbated SWI              
                C0i = Fnonbioi/((1-bsd.por)*bsd.w); % this value is not used anywhere - C0i is calculated from SWI-flux in benthic_zTOC_RCM.m
%                Fnonbioi = F.* ( swi.C0_nonbio*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
%                C0i = F.*swi.C0_nonbio;

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

        
    end
    
end
