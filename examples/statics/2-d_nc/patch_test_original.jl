using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked

println("Q4. Plane stress.")

E = 1.0
nu = 1.0 / 3
alpha, beta, gamma, delta, eta, phi =
    1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
ux(x, y) = alpha + beta * x + gamma * y
uy(x, y) = delta + eta * x + phi * y
MR = DeforModelRed2DStress

fens = FENodeSet(
    [
        1.0 -0.3
        2.3 -0.3
        2.3 0.95
        1.0 0.95
        1.4 0.05
        1.9 -0.03
        1.7 0.5
        1.3 0.6
    ],
)
fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

# Apply prescribed displacements to exterior nodes
for i = 1:4
    setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
    setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
end

applyebc!(u)
numberdofs!(u)

material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

F = zeros(nfreedofs(u))
# F = nzebcloadsstiffness(femm, geom, u)
K = stiffness(femm, geom, u)
K_ff, K_fd = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[(:ff, :fd)]
F_f = vector_blocked(F, nfreedofs(u))[:f]
U_d = gathersysvec(u, :d)
F_f = F_f - K_fd * U_d
K_ff = cholesky(K_ff)
U = K_ff \ (F_f)
scattersysvec!(u, U[:])

for i = 5:8
    uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
    println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
end

File = "a.vtk"
vtkexportmesh(
    File,
    fes.conn,
    geom.values,
    FinEtools.MeshExportModule.VTK.Q4;
    vectors = [("u", u.values)],
)

true

