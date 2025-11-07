using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
include("utilities.jl")


println("Q4. Plane stress.")

E = 1.0e5
nu = 1/3
MR = DeforModelRed2DStress
material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

width1 = 1.0
height1 = 2.0
n1x = 1
n1y = 2
fens1, fes1 = T3block(width1, height1, n1x, n1y)

boundaryfes1 = meshboundary(fes1)
edge_fes1 = subset(boundaryfes1, selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1], inflate=1e-8))


geom1 = NodalField(fens1.xyz)
u1 = NodalField(zeros(size(fens1.xyz, 1), 2)) # displacement field

box1 = [0.0,0.0,0.0,height1]
dbc_nodes1 = selectnode(fens1; box=box1, inflate=1e-8)
for i in dbc_nodes1
    setebc!(u1, [i], 1, 0.0)
    setebc!(u1, [i], 2, 0.0)
end

applyebc!(u1)
numberdofs!(u1)
femm = FEMMDeforLinear(MR, IntegDomain(fes1, GaussRule(2, 2)), material)

K1 = stiffness(femm, geom1, u1)
K1_ff = matrix_blocked(K1, nfreedofs(u1), nfreedofs(u1))[:ff]
F1 = zeros(size(K1, 1))
F1_ff = vector_blocked(F1, nfreedofs(u1))[:f]

edge_nodes1 = selectnode(fens1; box=[width1,width1, 0.0,height1], inflate=1e-8)
#########################################################################################
width2 = 1.0
height2 = 1.0
n2x = 1
n2y = 3
fens2, fes2 = T3block(width2, height2, n2x, n2y)
# shift the second mesh to the right by 1.0
fens2.xyz[:, 1] .+= 1.0

boundaryfes2 = meshboundary(fes2)
edge_fes2 = subset(boundaryfes2, selectelem(fens2, boundaryfes2, box=[1.0,1.0, 0.0,height2], inflate=1e-8))

geom2 = NodalField(fens2.xyz)
u2 = NodalField(zeros(size(fens2.xyz, 1), 2)) # displacement field
numberdofs!(u2)
femm2 = FEMMDeforLinear(MR, IntegDomain(fes2, GaussRule(2, 2)), material)
K2 = stiffness(femm2, geom2, u2)
K2_ff = matrix_blocked(K2, nfreedofs(u2), nfreedofs(u2))[:ff]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(u2))[:f]
edge_nodes2 = selectnode(fens2; box=[1.0,1.0, 0.0,height2], inflate=1e-8)
##########################################################################################

xs_i = [1.0,1.0,1.0]
ys_i = [0.0,1.0, 2.0]
fens_i, fes_i = L2blockx2D(xs_i, ys_i)
geom_i = NodalField(fens_i.xyz)
u_i = NodalField(zeros(size(fens_i.xyz, 1), 2)) # displacement field
numberdofs!(u_i)
D1,Pi_NC1,Pi_phi1 = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=1,tol=1e-8)
D2,Pi_NC2,Pi_phi2 = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=1,tol=1e-8)

