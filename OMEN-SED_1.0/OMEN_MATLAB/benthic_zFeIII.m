classdef benthic_zFeIII
    % Solve FeIII
    
    properties
        qdispFeIII=0.0;     %157.68;            % FeIII diffusion coefficient in water (cm2/yr)
        adispFeIII=0.0;     %7.884;             % FeIII linear coefficient for temperature dependence (cm2/yr/oC)
        DFeIII1;                      % FeIII diffusion coefficient in bioturbated layer (cm2/yr)
        DFeIII2;                      % FeIII diffusion coefficient in non-bioturbated layer (cm2/yr)
        
        reac1;
        reac2;
    end
    
    methods
        function obj = benthic_zFeIII(bsd, swi)
            obj.DFeIII1=bsd.Dbio;  	% FeIII diffusion coefficient in bioturbated layer (cm2/yr)
            obj.DFeIII2=0.009;  %(obj.qdispFeIII+obj.adispFeIII*swi.T).*bsd.dispFactor;  %      % FeIII diffusion coefficient in non-bioturbated layer (cm2/yr) - set as very small value so I can use the GBCM to solve it
            
            %reactive terms: OM degradation
            obj.reac1=-bsd.FeIIIC;
            obj.reac2=-bsd.FeIIIC;
            
        end
        
        function [flxswiFeIII] = calc_input(obj, bsd, swi, r, Fe3_conc_in)
            % Iteratively solve for correct SWI-concentration of FeIII+ (i.e. so influx is as specified)
            
            swi.FeIII0 = Fe3_conc_in;
            
            result = obj.calc(bsd, swi, r);
            flxswiFeIII = result.flxswiFeIII;
           
        end
        
        function r = calc(obj, bsd, swi, r)
                       
            % Iteratively solve for zfeIII
            
            % try zero flux at zinf and see if we have any FeIII left, also
            % calculate [FeIII] at zinf for advective loss
            [flxzfeIII, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r, 2);
            
            
            if r.zno3 == bsd.zinf
                r.zfeIII = bsd.zinf;
                bctype = 2;
            else
                
                fun=@(zfeIII) obj.calcbc(zfeIII,bsd,swi,r,1);   % FeIII consumption from Fe2 from below zfeIII:  - obj.calcFFeIII(zfeIII,bsd, swi, r);
%                fun=@(zfeIII) obj.calcbc(zfeIII,bsd,swi,r,1) + obj.calcFFeIII(zfeIII,bsd, swi, r);   % Here, include FeIII consumption from oxidation of H2S from below zfeIII:  - obj.calcFFeIII(zfeIII,bsd, swi, r);

                %             % try zero flux at zinf and see if we have any FeIII left
                %             [flxzfeIII, conczinf, flxswi,rtmp] = obj.calcbc(bsd.zinf, bsd, swi, r, 2);
                if bsd.usescalarcode
                    if conczinf >=0
                        r.zfeIII = bsd.zinf;
                        bctype = 2;
                    else
                        bctype = 1;
                        conczinf = 0.0;
                        [flxzfeIII, conczno3, flxswi,rtmp] = obj.calcbc(r.zno3, bsd, swi, r, 1);
                        %                        r.zfeIII=fzero(fun,[max(r.zno3, 1e-10), bsd.zbio],bsd.fzerooptions);    % Fe- reduction only in bioturbated zone assume Fe-reduction occurs just in bioturbated layer
                        r.zfeIII=fzero(fun,[max(r.zno3, 1e-10), bsd.zinf],bsd.fzerooptions);    % use this when diffusion <> 0, without diffusion zinf does not work as fewer int. const. are used in the lowest layer
                    end
                else  % vectorized version
                    bctype = (conczinf < 0)*1 + (conczfeIII>=0)*2;
                    zfeIII=fzero_vec(fun,max(r.zno3, 1e-10), bsd.zinf,bsd.fzerooptions);
                    r.zfeIII = (bctype==1).*zfeIII + (bctype==2).*bsd.zinf;
                end
                
            end
            [flxzfeIII, conczfeIII, flxswiFeIII, r] = obj.calcbc(r.zfeIII, bsd, swi, r, bctype);    % Dom18.05.2016: not necessary for bctype 2 (done in line 32 already)
            
            flxswiFeIII = flxswiFeIII - (1-bsd.por).*bsd.w.*(swi.FeIII0-conczinf);
            if(abs(flxswiFeIII) <= bsd.tol_const)
                flxswiFeIII = 0.0;
            end
            
            r.flxzfeIII = flxzfeIII;
            r.conczfeIII = conczfeIII;
            r.flxswiFeIII = flxswiFeIII;
        end
        
        function [flxzfeIII, conczfeIII, flxswi,r] = calcbc(obj, zfeIII, bsd, swi, r, bctype)
            % Calculate trial solution for given zfeIII, matching boundary conditions from layer-by-layer solutions
            
            if(swi.TwoG_OM_model)
                % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
                % layer 1: 0 < z < zox, passive diffn
                %      ls =      prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
                rFeIII.ls1 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,         0,     0, r.zox, obj.DFeIII1, obj.DFeIII2);
                % layer 2: zox < z < zno3, passive diffn
                rFeIII.ls2 = r.zTOC.prepfg_l12(bsd, swi, r, 0,         0,         0,  r.zox, r.zno3, obj.DFeIII1, obj.DFeIII2);
                % layer 3: zno3 < z < zfeIII, FeIII consumption by OM oxidation
                % DH: 15.07.20 guess the matching at zbio here has to be updated -- as 1 and 2 IC are matched!
                rFeIII.ls3 = r.zTOC.prepfg_l12(bsd, swi, r, obj.reac1, obj.reac2, 0, r.zno3, zfeIII, obj.DFeIII1, obj.DFeIII2);
                
                % Work up from the bottom, matching solutions at boundaries
                % Basis functions at bottom of layer 3 zfeIII
                [ e3_zfeIII, dedz3_zfeIII, f3_zfeIII, dfdz3_zfeIII, g3_zfeIII, dgdz3_zfeIII] ...
                    = r.zTOC.calcfg_l12(zfeIII, bsd, swi, r, obj.reac1, obj.reac2, 0, rFeIII.ls3);
                
                % Match at zno3, layer 2 - layer 3 (continuity and flux)
                % basis functions at bottom of layer 2
                [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                    = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r,     0,            0, 0, rFeIII.ls2);
                % ... and top of layer 3
                [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                    = r.zTOC.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1, obj.reac2, 0, rFeIII.ls3);
                % match solutions at zno3 - continuous concentration and flux
                [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                    e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                    0, 0);
                
                %%%%%%%%% DH 07.2020: TODO adjust here as well
                
                % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
                % flux of Fe2 to oxic interface (Source of FeIII)
                %            FH2S = 0.0; %r.zTOC.calcReac(r.zno3, zfeIII, bsd.FeIIIC, bsd.FeIIIC, bsd, swi, r) + 0.0; % no secondary redox!
                FFe2 = r.zTOC.calcReac(r.zno3, zfeIII, bsd.FeIIIC, bsd.FeIIIC, bsd, swi, r);
                %                + bsd.gammaCH4.*r.zTOC.calcReac(zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
                % basis functions at bottom of layer 1
                [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                    = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0 , 0 , 0, rFeIII.ls1);
                % basis functions at top of layer 2
                [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                    = r.zTOC.calcfg_l12(r.zox, bsd, swi, r, 0, 0, 0, rFeIII.ls2);
                % transform to use coeffs from l3
                [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                    zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
                
                % match solutions at zox - continuous concentration, flux discontinuity from H2S ox
                D = (r.zox <= bsd.zbio).*obj.DFeIII1 + (r.zox > bsd.zbio).*obj.DFeIII2;
                
                [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                    e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                    0, -r.zxf.*bsd.gammaFe2*FFe2./D);
                % with pyrite:      0, -r.zxf.*bsd.gammaFe2*FFe2./D);
               	% without pyrite:   0, -r.zxf.*bsd.gammaFe2.*FFe2./D);
                
                % Solution at swi, top of layer 1
                [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                    = r.zTOC.calcfg_l12(0, bsd, swi, r, 0 , 0 , 0, rFeIII.ls1);
                % transform to use coeffs from l3
                [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                    zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            else
                % Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
                % layer 1: 0 < z < zox, passive diffn
                rFeIII.ls1 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, 0,               0,     0, r.zox, obj.DFeIII1, obj.DFeIII2);
                % layer 2: zox < z < zno3, passive diffn
                rFeIII.ls2 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, 0,               0,  r.zox, r.zno3, obj.DFeIII1, obj.DFeIII2);
                % layer 3: zno3 < z < zfeIII, FeIII consumption by OM oxidation
                % DH: 15.07.20 guess the matching at zbio here has to be updated -- as 1 and 2 IC are matched!
                rFeIII.ls3 = r.zTOC_RCM.prepfg_l12(bsd, swi, r, obj.reac1, 0, r.zno3, zfeIII, obj.DFeIII1, obj.DFeIII2);
                
                % Work up from the bottom, matching solutions at boundaries
                % Basis functions at bottom of layer 3 zfeIII
                [ e3_zfeIII, dedz3_zfeIII, f3_zfeIII, dfdz3_zfeIII, g3_zfeIII, dgdz3_zfeIII] ...
                    = r.zTOC_RCM.calcfg_l12(zfeIII, bsd, swi, r, obj.reac1, 0, rFeIII.ls3);
%                 if(isinf(f3_zfeIII) || isinf(dfdz3_zfeIII))
%                     f3_zfeIII = realmax;
%                     dfdz3_zfeIII = realmax;
%                 end
                
                % Match at zno3, layer 2 - layer 3 (continuity and flux)
                % basis functions at bottom of layer 2
                [ e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3] ...
                    = r.zTOC_RCM.calcfg_l12(r.zno3, bsd, swi, r,     0, 0, rFeIII.ls2);
                % ... and top of layer 3
                [ e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3] ...
                    = r.zTOC_RCM.calcfg_l12(r.zno3, bsd, swi, r, obj.reac1, 0, rFeIII.ls3);
                % match solutions at zno3 - continuous concentration and flux
                [zno3.a, zno3.b, zno3.c, zno3.d, zno3.e, zno3.f] = benthic_utils.matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, ...
                    e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, ...
                    0, 0);
                
                %%%%%%%%% DH 07.2020: TODO adjust here as well
                
                % Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
                % flux of Fe2 to oxic interface (Source of FeIII)
                %            FH2S = 0.0; %r.zTOC_RCM.calcReac(r.zno3, zfeIII, bsd.FeIIIC, bsd.FeIIIC, bsd, swi, r) + 0.0; % no secondary redox!
                FFe2 = r.zTOC_RCM.calcReac(r.zno3, zfeIII, bsd.FeIIIC, bsd, swi, r);
                %                + bsd.gammaCH4.*r.zTOC_RCM.calcReac(zfeIII, bsd.zinf, bsd.MC, bsd.MC, bsd, swi, r); % Dominik 25.02.2016
                % basis functions at bottom of layer 1
                [ e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox] ...
                    = r.zTOC_RCM.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rFeIII.ls1);
                % basis functions at top of layer 2
                [ e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox] ...
                    = r.zTOC_RCM.calcfg_l12(r.zox, bsd, swi, r, 0, 0, rFeIII.ls2);
                % transform to use coeffs from l3
                [e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox] = benthic_utils.xformsoln(e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                    zno3.a , zno3.b , zno3.c , zno3.d , zno3.e ,zno3.f);
                
                % match solutions at zox - continuous concentration, flux discontinuity from H2S ox
                D = (r.zox <= bsd.zbio).*obj.DFeIII1 + (r.zox > bsd.zbio).*obj.DFeIII2;
                
                [zox.a, zox.b, zox.c, zox.d, zox.e, zox.f] = benthic_utils.matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, ...
                    e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, ...
                    0, -r.zxf.*bsd.gammaFe2*FFe2./D);
                % with pyrite:      0, -r.zxf.*bsd.gammaFe2*FFe2./D);
               	% without pyrite:   0, -r.zxf.*bsd.gammaFe2.*FFe2./D);                

                % Solution at swi, top of layer 1
                [ e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0] ...
                    = r.zTOC_RCM.calcfg_l12(0, bsd, swi, r, 0 , 0, rFeIII.ls1);
                % transform to use coeffs from l3
                [ e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0]= benthic_utils.xformsoln(e1_0, f1_0, g1_0, dedz1_0, dfdz1_0, dgdz1_0, ...
                    zox.a , zox.b , zox.c , zox.d , zox.e ,zox.f);
            end
            %%%%%%%%%  DH 07.2020: TODO here equations must solved differently, as
            %%%%%%%%%  viewer BC (compare PO_M)
            
            % Find solutions for two possible types of lower bc
            %  case 1  zero concentration at zfeIII
            % Solve for AFeIII, BFeIII given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
            % AFeIII*e3_zfeIII   +  BFeIII*f3_zfeIII  + g3_zfeIII = 0;
            % AFeIII*e1_0     +   BFeIII*f1_0     + g1_0  = swi.FeIII0;
            
            % | e3_zfeIII f3_zfeIII |  |AFeIII|   = | -g3_zfeIII       |
            % | e1_0     f1_0   |  |BFeIII|     | swi.FeIII0 - g1_0 |
            [ bctype1_A3, bctype1_B3]      = benthic_utils.solve2eqn(e3_zfeIII, f3_zfeIII, e1_0, f1_0, -g3_zfeIII, swi.FeIII0 - g1_0);
            
            A3_test1 = -g3_zfeIII/e3_zfeIII;
            A3_test2 = (swi.FeIII0 - g1_0)/e1_0;
            
            %  case  2 zero flux at zfeIII
            % AFeIII*de3dz_zfeIII   +  BFeIII*dfdz3_zfeIII  + dgdz3_zfeIII = 0;
            % AFeIII*e1_0         +   BFeIII*f1_0       + g1_0       = swi.FeIII0;
            [ bctype2_A3, bctype2_B3]      = benthic_utils.solve2eqn(dedz3_zfeIII, dfdz3_zfeIII, e1_0, f1_0, -dgdz3_zfeIII, swi.FeIII0 - g1_0);
            
            % Choose type of solution requested (vectorized form)
            rFeIII.A3 = (bctype==1).*bctype1_A3 + (bctype==2).*bctype2_A3;
            rFeIII.B3 = (bctype==1).*bctype1_B3 + (bctype==2).*bctype2_B3;
            
            %             rFeIII.A3 = A3_test1;   %(bctype==1).*bctype1_A3 + (bctype==2).*bctype2_A3;
            %             rFeIII.B3 = 0.0;        % (bctype==1).*bctype1_B3 + (bctype==2).*bctype2_B3;
            
            % calculate conc and flux at zfeIII
            conczfeIII = rFeIII.A3.*e3_zfeIII+rFeIII.B3.*f3_zfeIII + g3_zfeIII;
            D = (zfeIII <= bsd.zbio).*obj.DFeIII1 + (zfeIII > bsd.zbio).*obj.DFeIII2;
            flxzfeIII = D.*(rFeIII.A3.*dedz3_zfeIII+rFeIII.B3.*dfdz3_zfeIII + dgdz3_zfeIII);        % includes 1/por ie flux per (cm^2 pore area)
            
            % flux at swi - DO include por so this is per cm^2 water column area
            % DH: added advective flux 28.05.2016
            flxswi = (1-bsd.por).*(obj.DFeIII1.*(rFeIII.A3.*dedz1_0+rFeIII.B3.*dfdz1_0 + dgdz1_0)); % - bsd.w.*swi.FeIII0);   % NB: use A3, B3 as these are _xformed_ layer 1 basis functions
            
            % save coeffs for layers 2 and 1
            rFeIII.A2 = zno3.a.*rFeIII.A3 + zno3.b.*rFeIII.B3 + zno3.e;
            rFeIII.B2 = zno3.c.*rFeIII.A3 + zno3.d.*rFeIII.B3 + zno3.f;
            
            rFeIII.A1 = zox.a.*rFeIII.A3 + zox.b.*rFeIII.B3 + zox.e;
            rFeIII.B1 = zox.c.*rFeIII.A3 + zox.d.*rFeIII.B3 + zox.f;
            
            
            r.rFeIII = rFeIII;
            
        end
        
        
        function FFeIII = calcFFeIII(obj, zfeIII, bsd, swi, r)
            %% Calculate FeIII consumption below zfeIII, by organic matter and indirectly via methane oxidation
            % no need to simulate FeIII consumption from Fe2 from below zfeIII (bc it is not produced there)
            % but FeIII is used to oxidize some of the H2S produced!
            tmpreac1    = bsd.gammaH2SFe*bsd.FeIIIH2S*bsd.SO4C;
            tmpreac2    = bsd.gammaH2SFe*bsd.FeIIIH2S*bsd.SO4C;
            %             FFeIII = r.zTOC.calcReac(zfeIII, bsd.zinf, tmpreac1, tmpreac2, bsd, swi, r);
            
            if(swi.TwoG_OM_model)
                FFeIII = r.zTOC.calcReac(zfeIII, bsd.zinf, tmpreac1, tmpreac2, bsd, swi, r);
            else
                FFeIII = r.zTOC_RCM.calcReac(zfeIII, bsd.zinf, tmpreac1, bsd, swi, r);
            end
            FFeIII = 0.0;    % no oxidation of H2S by FeIII
        end
        

        
        function [FeIII, flxFeIII] = calcFeIII(obj, z, bsd, swi, r)
            % Calculate FeIII concentration and flux at depth z from solution
            
            rFeIII = r.rFeIII;
            if z <= r.zfeIII
                if z <= bsd.zbio
                    D = obj.DFeIII1;
                else
                    D = obj.DFeIII2;
                end
                
                if(swi.TwoG_OM_model)
                    if z <= r.zox   % layer 1
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFeIII.ls1);
                        FeIII     = r.rFeIII.A1.*e + r.rFeIII.B1.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A1.*dedz+r.rFeIII.B1.*dfdz + dgdz);
                    elseif z <= r.zno3 % layer 2
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFeIII.ls2);
                        FeIII     = r.rFeIII.A2.*e + r.rFeIII.B2.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A2.*dedz+r.rFeIII.B2.*dfdz + dgdz);
                    else
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1, obj.reac2 , 0, rFeIII.ls3);
                        FeIII     = r.rFeIII.A3.*e + r.rFeIII.B3.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A3.*dedz+r.rFeIII.B3.*dfdz + dgdz);
                    end
                    
                else
                    if z <= r.zox   % layer 1
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0, 0, rFeIII.ls1);
                        FeIII     = r.rFeIII.A1.*e + r.rFeIII.B1.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A1.*dedz+r.rFeIII.B1.*dfdz + dgdz);
                    elseif z <= r.zno3 % layer 2
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0 , 0, rFeIII.ls2);
                        FeIII     = r.rFeIII.A2.*e + r.rFeIII.B2.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A2.*dedz+r.rFeIII.B2.*dfdz + dgdz);
                    else
                        [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, obj.reac1, 0, rFeIII.ls3);
                        FeIII     = r.rFeIII.A3.*e + r.rFeIII.B3.*f + g;
                        flxFeIII  = D.*(r.rFeIII.A3.*dedz+r.rFeIII.B3.*dfdz + dgdz);
                    end
                    
                end
                
            else
                
                FeIII = 0;
                flxFeIII = 0;
            end
        end
        
    	function [Fe3_con] = calcConcentration(obj, bsd, swi, r, Fe3_flux)
            % Calculate FeIII SWI-concentration (mol/cm^3) from flux (mol/(cm^2 yr)) at SWI for the case zno3 < zbio
            
            % First - BC 4: calculate g & dgdz (i.e. PhiI, PhiII, PhiIII) coming from POC degradation from Fe-reduction at zbio            
            [ e_zbio, dedz_zbio, f_zbio, dfdz_zbio, g_zbio, dgdz_zbio]  = r.zTOC_RCM.calcfg_l1( bsd.zbio, bsd, swi, r, obj.reac1, obj.DFeIII1, 0.0);
            
            % Second - BC 3: calculate g & dgdz (i.e. PhiI, PhiII, PhiIII) coming from POC degradation from Fe-reduction at zno3            
            [ e_zno3, dedz_zno3, f_zno3, dfdz_zno3, g_zno3, dgdz_zno3]  = r.zTOC_RCM.calcfg_l1( r.zno3, bsd, swi, r, obj.reac1, obj.DFeIII1, 0.0);
            
            b1 = bsd.w/obj.DFeIII1;
            b2 = bsd.w/obj.DFeIII1;
            
            % solution obtained using Sage by solving the Flux divergence equation (see Eq. (51 in GMD paper) for Fe3+ by using BC3 and BC4
            Fe3_con =  (obj.DFeIII1*dgdz_zno3*bsd.por -  obj.DFeIII1*dgdz_zno3 - (obj.DFeIII1*dgdz_zbio*bsd.por - obj.DFeIII1*dgdz_zbio)*exp(-b2*bsd.zbio + b2*r.zno3) - Fe3_flux*exp(b1*r.zno3))*exp(-b1*r.zno3)/((bsd.por - 1)*bsd.w);

%             rFeIII = r.rFeIII;
%             if z <= r.zfeIII
%                 if z <= bsd.zbio
%                     D = obj.DFeIII1;
%                 else
%                     D = obj.DFeIII2;
%                 end
%                 
%                 if(swi.TwoG_OM_model)
%                     if z <= r.zox   % layer 1
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFeIII.ls1);
%                         FeIII     = r.rFeIII.A1.*e + r.rFeIII.B1.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A1.*dedz+r.rFeIII.B1.*dfdz + dgdz);
%                     elseif z <= r.zno3 % layer 2
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, 0 , 0 , 0, rFeIII.ls2);
%                         FeIII     = r.rFeIII.A2.*e + r.rFeIII.B2.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A2.*dedz+r.rFeIII.B2.*dfdz + dgdz);
%                     else
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC.calcfg_l12(z, bsd, swi, r, obj.reac1, obj.reac2 , 0, rFeIII.ls3);
%                         FeIII     = r.rFeIII.A3.*e + r.rFeIII.B3.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A3.*dedz+r.rFeIII.B3.*dfdz + dgdz);
%                     end
%                     
%                 else
%                     if z <= r.zox   % layer 1
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0, 0, rFeIII.ls1);
%                         FeIII     = r.rFeIII.A1.*e + r.rFeIII.B1.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A1.*dedz+r.rFeIII.B1.*dfdz + dgdz);
%                     elseif z <= r.zno3 % layer 2
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, 0 , 0, rFeIII.ls2);
%                         FeIII     = r.rFeIII.A2.*e + r.rFeIII.B2.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A2.*dedz+r.rFeIII.B2.*dfdz + dgdz);
%                     else
%                         [ e, dedz, f, dfdz, g, dgdz]  = r.zTOC_RCM.calcfg_l12(z, bsd, swi, r, obj.reac1, 0, rFeIII.ls3);
%                         FeIII     = r.rFeIII.A3.*e + r.rFeIII.B3.*f + g;
%                         flxFeIII  = D.*(r.rFeIII.A3.*dedz+r.rFeIII.B3.*dfdz + dgdz);
%                     end
%                     
%                 end
%                 
%             else
%                 
%                 FeIII = 0;
%                 flxFeIII = 0;
%             end
        end
        
    end
    
end

