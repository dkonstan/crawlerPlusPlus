import MDAnalysis as mda

# TIP3P
COULOMBS_PER_E = 1.60218e-19
KG_PER_AMU = 1.66054e-27
massO = 15.99943 * KG_PER_AMU
massH = 1.007947 * KG_PER_AMU
massLi = 6.941 * KG_PER_AMU
chargeO = -0.834 * COULOMBS_PER_E
chargeH = 0.417 * COULOMBS_PER_E
chargeLi = 1.0 * COULOMBS_PER_E
rOH = 0.9572e-10
kOH = 462750.4 * (1000) * (1 / 6.0221408e23) * (1e9)**2  # originally in kJ/mol/nm^2
kAngle = 836.8 * (1000) * (1 / 6.0221408e23)  # originally in kJ/mol/radian^2
liEpsilon = 0.1171084864 * (1000) * (1 / 6.0221408e23)  # originally in kJ/mol
liSigma = 0.18263423721876956e-9  # originally in nm
oEpsilon = 0.635968 * (1000) * (1 / 6.0221408e23)  # originally inkJ/mol
oSigma = 0.31507524065751241e-9  # originally in nm
hEpsilon = 0.0 * (1000) * (1 / 6.0221408e23)  # originally in kJ/mol
hSigma = 1.0e-9  # originally in nm

system = mda.Universe("Li_waterClusterReordered.xyz")
waters = system.select_atoms("name O H")

with open("Li+W20_tip3p.top", "w+") as top:

    top.write("<atomtypes>\n")
    for idx, atom in enumerate(system.atoms):
        top.write(f"{idx} {atom.name}\n")
    top.write("<end>\n\n")

    top.write("<masses>\n")
    cnt = 0
    for i in range(waters.n_atoms // 3):
        top.write(f"{cnt} {massO}\n")
        top.write(f"{cnt + 1} {massH}\n")
        top.write(f"{cnt + 2} {massH}\n")
        cnt += 3
    top.write(f"{cnt} {massLi}\n")
    top.write("<end>\n\n")

    top.write("<charges>\n")
    cnt = 0
    for i in range(waters.n_atoms // 3):
        top.write(f"{cnt} {chargeO}\n")
        top.write(f"{cnt + 1} {chargeH}\n")
        top.write(f"{cnt + 2} {chargeH}\n")
        cnt += 3
    top.write(f"{cnt} {chargeLi}\n")
    top.write("<end>\n\n")

    top.write("<bonds>\n")
    for i in range(0, waters.n_atoms, 3):
        top.write(f"{i} {i + 1} {rOH} {kOH}\n")
        top.write(f"{i} {i + 2} {rOH} {kOH}\n")
    top.write("<end>\n\n")

    top.write("<angles>\n")
    for i in range(0, waters.n_atoms, 3):
        top.write(f"{i + 1} {i} {i + 2} 104.5 {kAngle}\n")
    top.write("<end>\n\n")

    top.write("<vdw>\n")
    cnt = 0
    for i in range(0, waters.n_atoms, 3):
        top.write(f"{i} {oSigma} {oEpsilon}\n")
        top.write(f"{i + 1} {hSigma} {hEpsilon}\n")
        top.write(f"{i + 2} {hSigma} {hEpsilon}\n")
        cnt += 3
    top.write(f"{cnt} {liSigma} {liEpsilon}\n")
    top.write("<end>\n\n")

    top.write("<nonbonded_exclude>\n")
    for i in range(0, waters.n_atoms, 3):
        top.write(f"{i} {i + 1}\n")
        top.write(f"{i} {i + 2}\n")
        top.write(f"{i + 1} {i + 2}\n")
    top.write("<end>\n\n")

    top.write("<box>\n")
    top.write("1.0 1.0 1.0\n")
    top.write("<end>\n\n")
