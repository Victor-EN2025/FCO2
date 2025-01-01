function [F_CO2, dpCO2, K_CO2] = FCO2v2(pCO2_water, pCO2_atm, SST, SSS, u)
   %%Function for the calculation or air-sea CO2 fluxes
        % Victor Ebolo Nkongo
        %University of Douala/ Instute of Fisheries and Aquatics Sciences 
        %   filename: FCO2v2.m
        % Current:
        %   FCO2V2 ( 07/06/2024 by Victor Ebolo Nkongo)
        %   FCO2 ( 03/23/2015 by Cecilia Chapa-Balcorta)
        %   copyright @ 2024 Victor Ebolo Nkongo 
       
        %INPUT
    
    % Inputs:
    % pCO2_water - seawater pCO2 (uatm)
    % pCO2_atm - atmospheric pCO2 (uatm)
    % SST - Temperature (Celsius)
    % S - Salinity
    % u - Wind speed (m/s)
    
    %Syntaxe: 
    %[Fluxe,dpco2,Solubility]=FCO2V2(pCO2_water,xCO2_atm,SST,SSS,wing);

    % Calculate the Schmidt number for CO2 in seawater
    Sc = Schmidt(SST);

    % Calculate the gas transfer velocity K using the updated coefficient
    K_CO2 = 0.251 .* (u.^2) .* ((Sc / 660) .^ -0.5);

    % Calculate the difference in pCO2 between water and atmosphere
    dpCO2 = pCO2_water - pCO2_atm;

    % Calculate CO2 solubility using the Weiss (1974) formula
    a = Ko_weiss(SST, SSS);

    % Calculate the air-sea CO2 flux in mmol m^-2 d^-1
    F_CO2 = 0.24 * K_CO2 .* a .* dpCO2;
end
 %%%%%Subrutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%****Solibuility constant (Weiss, 1974) ***********************************

function [Ko]=Ko_weiss(SST,SSS)
%
A=[-58.0931, 90.5069, 22.2940];  %mol/Kg.atm
B=[0.027766, -0.025888, 0.0050578]; %mol/Kg.atm
SST=SST+273.15; %Conversio from Celsius degrees to Kelvins
Ln_Ko=A(1)+(A(2).*(100./SST))+(A(3).*log(SST./100))+SSS.*(B(1)+(B(2).*(SST./100))+(B(3).*(SST./100).^2));
Ko=exp(Ln_Ko);

end
     
%******** Schmidt Number*********

    %For water of salinity=35 and temperature range –2° and 40°C %%%%%%%%%%%%%
    
    function [Sc]=Schmidt(SST)
            A = 2116.8;     B = -136.25;     C = 4.7353;     D = -0.092307;
            E=0.0007555;
            Sc= A + (B.*SST)+(C.*SST.^2)+(D.*SST.^3) + (E.*SST.^4);
               
    end
        
