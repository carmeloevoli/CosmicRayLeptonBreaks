def pwn_paradigm():
    cLightCGS = 2.99792458e10 # cm / s 
    yearCGS = 3.1e7 # s
    kpcCGS = 3.11e21 # s

    E2I = 0.2 / 1e4 # GeV cm-2 s-1 sr-1
    epsilon_pos = 4. * math.pi / cLightCGS * E2I # GeV / cm3
    print (f'energy density positrons : {epsilon_pos * 1e9:4.0e} eV/cm3')

    epsilon_pos /= GeV2erg # erg / cm3
    print (f'energy density positrons : {epsilon_pos:4.0e} erg/cm3')

    tau_loss = 4. * 1e6 * yearCGS
    R_G = 15. * kpcCGS
    l_G = 5. * kpcCGS
    V_G = 2. * math.pi * R_G**2.0 * l_G

    print(f'V_G : {np.power(V_G, 1. / 3.) / kpcCGS:4.0f} kpc3')

    print(f'L pos : {epsilon_pos * V_G / tau_loss:4.0e} erg/s')

    E_PWN = 5e47 # erg
    R_PWN = 1. / 50. / yearCGS

    print(f'L PWN : {E_PWN * R_PWN:4.0e} erg/s')