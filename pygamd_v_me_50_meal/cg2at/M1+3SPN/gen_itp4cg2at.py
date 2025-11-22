rtp = """P       P2  1.5000   1
      O1P      ON3 -0.7800   2
      O2P      ON3 -0.7800   2
      C5'     CN8B -0.0800   2
      H5'      HN8  0.0900   2
     H5''      HN8  0.0900   2
      O5'      ON2 -0.5700   2
      C4'      CN7  0.1600   3
      H4'      HN7  0.0900   3
      O4'      ON6 -0.5000   3
      C3'      CN7  0.0100   6
      H3'      HN7  0.0900   6
      O3'      ON2 -0.5700   6
      C2'      CN8 -0.1800   5
      H2'      HN8  0.0900   5
     H2''      HN8  0.0900   5
      C1'     CN7B  0.1600   3
      H1'      HN7  0.0900   3
       N1     NN2B -0.3400   4
       C2     CN1T  0.5100   4
       O2      ON1 -0.4100   4
       N3     NN2U -0.4600   4
       H3      HN2  0.3600   4
       C4      CN1  0.5000   4
       O4      ON1 -0.4500   4
       C5     CN3T -0.1500   4
      C5M      CN9 -0.1100   4
      H51      HN9  0.0700   4
      H52      HN9  0.0700   4
      H53      HN9  0.0700   4
       C6      CN3  0.1700   4
       H6      HN3  0.1700   4"""

bond = """P   O1P
        P   O2P
        P   O5'
      O5'   C5'
      C5'   C4'
      C4'   O4'
      C4'   C3'
      O4'   C1'
      C1'    N1
      C1'   C2'
       N1    C2
       N1    C6
       C2    N3
       N3    H3
       N3    C4
       C4    C5
       C5   C5M
      C2'   C3'
      C3'   O3'
      C1'   H1'
      C2'  H2''
      C3'   H3'
      C4'   H4'
      C5'   H5'
      C5'  H5''
       C6    H6
      C5M   H51
      C5M   H52
      C5M   H53
       C2    O2
       C4    O4
       C5    C6
      C2'   H2'"""
mass_dict_atom = {"C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, 'H': 1.008, 'S': 32}
q_tot = 0
lines = rtp.split('\n')
atoms = []
for line_idx in range(len(lines)):
    atom_name, ff_name, charge, _ = lines[line_idx].split()
    # print(atom_name, ff_name, charge)
    # if 'H' in atom_name:
    #     continue
    atoms.append(atom_name)
    q_tot += float(charge)
    # 1        CN7      1    ATP    C4'      1       0.16     13.019   ; qtot 0.16
    print(f"{line_idx+1:>6d} {ff_name:>10s}       1    Gb  {atom_name:>5s}      1     {float(charge):>6.2f}   {mass_dict_atom[atom_name[0]]:>8.3f}   ; qtot {q_tot:6.2f}")

bonds = []
for line in bond.split('\n'):
    atom1, atom2 = line.split()
    # if 'H' in atom1 or 'H' in atom2:
    #     continue
    cur_bond = (atoms.index(atom1), atoms.index(atom2))
    bonds.append(cur_bond)
    print(f"{cur_bond[0]+1:>6d} {cur_bond[1]+1:>6d}  1")
