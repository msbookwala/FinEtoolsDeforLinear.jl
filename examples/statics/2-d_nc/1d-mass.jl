using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked

E = 1.0e5
nu = 1/3
MR = DeforModelRed1DStrain

fens, fes = L2block(1.0, 3)
geom = NodalField(fens.xyz)
u = NodalField(zeros(count(fens), 2))
numberdofs!(u)
material = MatDeforElastIso(MR, 1.0, E, nu, 0.0)
femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(1, 2)), material)
M = mass(femm, geom, u)