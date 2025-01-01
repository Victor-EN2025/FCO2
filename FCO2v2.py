# === Fonction principale : Calcul des flux de CO2 ===
def FCO2v2(pCO2_water, pCO2_atm, SST, SSS, u):
    """
    Function for the calculation of air-sea CO2 fluxes
    Victor Ebolo Nkongo
    University of Douala/Institute of Fisheries and Aquatics Sciences
    
    filename: FCO2v2.py
    Current:
        FCO2v2 (Created on Mon Dec 30 12:33:21 2024 by Victor Ebolo Nkongo)
       
        copyright @ 2024 Victor Ebolo Nkongo
    
    Inputs:
        pCO2_water - seawater pCO2 (µatm)
        pCO2_atm - atmospheric pCO2 (µatm)
        SST - Temperature (Celsius)
        SSS - Salinity
        u - Wind speed (m/s)
    
    Syntax:
        [Fluxe, dpCO2, Solubility] = FCO2V2(pCO2_water, pCO2_atm, SST, SSS, u)
    
    Outputs:
        F_CO2 - Air-sea CO2 flux (mmol/m²/day)
        dpCO2 - Difference in pCO2 (µatm)
        K_CO2 - Gas transfer velocity (m/day)
    """
    # Calculate the Schmidt number for CO2 in seawater
    Sc = Schmidt(SST)

    # Calculate the gas transfer velocity KCO2 (wanninkhof, 2014) using the updated coefficient
    K_CO2 = 0.251 * (u**2) * ((Sc / 660) ** -0.5)

    # Calculate the difference in pCO2 between water and atmosphere
    dpCO2 = pCO2_water - pCO2_atm

    # Calculate CO2 solubility using the Weiss (1974) formula
    a = Ko_weiss(SST, SSS)

    # Calculate the air-sea CO2 flux in mmol/m²/day
    F_CO2 = 0.24 * K_CO2 * a * dpCO2

    return F_CO2, dpCO2, K_CO2


# === Subroutines ===

# Solubility constant (Weiss, 1974)
import numpy as np

def Ko_weiss(SST, SSS):
    """
    Calculate CO2 solubility using the Weiss (1974) formula.
    """
    A = [-58.0931, 90.5069, 22.2940]  # mol/kg.atm
    B = [0.027766, -0.025888, 0.0050578]  # mol/kg.atm
    SST = SST + 273.15  # Conversion from Celsius to Kelvin
    Ln_Ko = (
        A[0] + (A[1] * (100. / SST)) + (A[2] * np.log(SST / 100.)) +
        SSS * (B[0] + (B[1] * (SST / 100.)) + (B[2] * (SST / 100.)**2))
    )
    return np.exp(Ln_Ko)

# Schmidt Number calculation
def Schmidt(SST):
    """
    Calculate the Schmidt number for CO2 in seawater.
    For water of salinity = 35 and temperature range –2° to 40°C.
    """
    A = 2116.8
    B = -136.25
    C = 4.7353
    D = -0.092307
    E = 0.0007555
    return A + (B * SST) + (C * SST**2) + (D * SST**3) + (E * SST**4)
