using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked

println("Q4. Plane stress.")

E = 1.0e5
nu = 1/3
alpha, beta, gamma, delta, eta, phi =
    1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
ux(x, y) = alpha + beta * x + gamma * y
uy(x, y) = delta + eta * x + phi * y
MR = DeforModelRed2DStress

width1 = 1.0
height1 = 2.0
n1x = 1
n1y = 2
fens1, fes1 = T3block(width1, height1, n1x, n1y)


geom = NodalField(fens1.xyz)
u = NodalField(zeros(size(fens1.xyz, 1), 2)) # displacement field

# Apply prescribed displacements to exterior nodes
box1 = [0.0,0.0,0.0,height1]
dbc_nodes1 = selectnode(fens1; box=box1, inflate=1e-8)
for i in dbc_nodes1
    setebc!(u, [i], 1, 0.0)
    setebc!(u, [i], 2, 0.0)
end

applyebc!(u)
numberdofs!(u)

material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

femm = FEMMDeforLinear(MR, IntegDomain(fes1, GaussRule(2, 2)), material)

# F = nzebcloadsstiffness(femm, geom, u)
K = stiffness(femm, geom, u)
K_ff = matrix_blocked(K, nfreedofs(u))[:ff]
K_ff = cholesky(K_ff)
F = [0.5,0,1,0,0.5,0]
U = K_ff \ (F)
scattersysvec!(u, U[:])

File = "a.vtk"
vtkexportmesh(
    File,
    fes1.conn,
    geom.values,
    FinEtools.MeshExportModule.VTK.T3;
    vectors = [( "u", u.values)],
)


