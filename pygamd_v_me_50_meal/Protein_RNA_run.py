#!/usr/bin/python
import sys
import math

from poetry import cu_gala as gala
from poetry import force_field_gala
from poetry import _options

filename = "particles.0102000000.xml"
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.0001
app = gala.Application(all_info, dt)

Temperature = 310.0  # (k)
concent = 50.0  # salt concentration

kappaD = 13.603 * math.sqrt(
    (50.0 / concent) * (Temperature / 300.0))  # related to salt concentration 13.603 -> 50mM NaCl
print("kappaD = ", kappaD)
real_epsilon = 0.26 * 4184.0  # (0.26 kcal/mol = 0.26**4184.0 J/mol)
R = 8.314472  # gas constant
enegy_reduce_unit = 1000.0
lenscale = 10.0  # angstrom
lsq = lenscale * lenscale  #
epsilon = real_epsilon / enegy_reduce_unit
print("epsilon = ", epsilon)

rcut = 4.0

neighbor_list = gala.NeighborList(all_info, rcut, 0.4)  # (,rcut,rbuffer)
neighbor_list.exclusion(["bond"])
neighbor_list.addExclusionsFromBonds()
neighbor_list.countExclusions()

ShortRangeEpsilon = 0.8368
debye_length = 0.794
ahdh = force_field_gala.WFDHForce(all_info, neighbor_list, rcut,
                  debye_length, "new.force_field")
app.add(ahdh)

bondforce_pro = gala.BondForceHarmonic(all_info)
bond_type = all_info.getBondInfo().getBondTypes()
print(bond_type)
for i in bond_type:
    bondforce_pro.setParams(i, 8033.28, 0.38)

app.add(bondforce_pro)


Temperature = 310.000  # k
T = Temperature * 8.3143 / 1000.0  # reduced unit

groupall = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, groupall)

#gp = gala.LangevinNVT(all_info, groupall, T, 12345)
#app.add(gp)

#bgroup = gala.ParticleSet(all_info, 'body')
#rigidbd = gala.LangevinNVTRigid(all_info, bgroup, T, 12345)
#app.add(rigidbd)

nbgroup = gala.ParticleSet(all_info, 'all')
bd = gala.LangevinNVT(all_info, nbgroup, T, 12345)
app.add(bd)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(500)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(2000)  # (period)
app.add(ZeroMomentum)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(int(1000))
app.add(DInfo)

# dcd = gala.DCDDump(all_info, 'wrap-particles', True)
# dcd.unwrap(True)
# dcd.setPeriod(int(1e5))
# app.add(dcd)

xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(int(1e6))  # (period)
xml.setOutputBond(True)
xml.setOutputVelocity(True)
xml.setOutputMass(True)
xml.setOutputCharge(True)
xml.setOutputBody(True)
app.add(xml)

# ready ro run

app.run(int(1e6))  # (How many steps to run)

app.setDt(0.001)
app.run(int(1e6))

app.setDt(0.02)
app.run(int(1e8))  # (How many steps to run)
# neighbor_list.printStats()
