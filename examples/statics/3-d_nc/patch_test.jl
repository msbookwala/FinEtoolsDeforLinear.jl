using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using Infiltrator
include("meshrefine.jl")


println("Q4. Plane stress.")

E = 1.0
nu = 1/3
MR = DeforModelRed3D
material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

N_elem1 = 5
N_elem2 = 3
N_elem_i = min(N_elem1, N_elem2)
left_m = "t"
right_m = "t"
bend = 0.0
lam_order = 1


alpha0, alpha1, alpha2, alpha3, beta0, beta1, beta2, beta3, gamma0, gamma1, gamma2, gamma3 =
    1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35, 1.0 / 29, -1.0 / 31, -1.0 / 27, -1.0 / 33, -1.0 / 28, 1.0 / 32
ux(x, y, z) = alpha0 + alpha1 * x + alpha2 * y + alpha3 * z
uy(x, y, z) = beta0 + beta1 * x + beta2 * y + beta3 * z
uz(x, y, z) = gamma0 + gamma1 * x + gamma2 * y + gamma3 * z

#########################################################################################
width1 = 0.5
height1 = 1.0
depth1 = 1.0
if left_m == "h"
    fens1, fes1 = H8block(width1, height1, depth1, floor(Int, N_elem1), N_elem1, N_elem1)
    Rule1 = GaussRule(3,2)
else
    fens1, fes1 = T4block(width1, height1, depth1, floor(Int, N_elem1), N_elem1, N_elem1)
    Rule1 = TetRule(4)
end

boundaryfes1 = meshboundary(fes1)
edge_fe_idx1 = selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1, 0.0,depth1], inflate=1e-8)
edge_fes1 = subset(boundaryfes1, edge_fe_idx1)
dbc_fes_idx1 = setdiff(1:count(boundaryfes1), edge_fe_idx1)


fens1.xyz[:, 1] .+= bend * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 0.5).^2#


geom1 = NodalField(fens1.xyz)
u1 = NodalField(zeros(size(fens1.xyz, 1), 3)) # displacement field

dbc_nodes1 = unique(stack(boundaryfes1.conn[dbc_fes_idx1]))
for i in dbc_nodes1
    setebc!(u1, [i], 1, ux(fens1.xyz[i, :]...))
    setebc!(u1, [i], 2, uy(fens1.xyz[i, :]...))
    setebc!(u1, [i], 3, uz(fens1.xyz[i, :]...))
end

applyebc!(u1)
numberdofs!(u1)
femm = FEMMDeforLinear(MR, IntegDomain(fes1, Rule1), material)

K1 = stiffness(femm, geom1, u1)
K1_ff = matrix_blocked(K1, nfreedofs(u1), nfreedofs(u1))[:ff]
K1_fd = matrix_blocked(K1, nfreedofs(u1), nfreedofs(u1))[:fd]
F1 = zeros(size(K1, 1))
F1_ff = vector_blocked(F1, nfreedofs(u1))[:f] - K1_fd * gathersysvec(u1, :d)

edge_nodes1 = selectnode(fens1; box=[width1,width1, 0.0,height1, 0.0,depth1], inflate=1e-8)
#########################################################################################
width2 = 0.5
height2 = 1.0
depth2 = 1.0
if right_m == "h"
    fens2, fes2 = H8block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = GaussRule(3,2)
else
    fens2, fes2 = T4block(width2, height2, depth2, floor(Int, N_elem2), N_elem2, N_elem2)
    Rule2 = TetRule(4)
end 
# shift the second mesh to the right by 1.0
fens2.xyz[:, 1] .+= 0.5


boundaryfes2 = meshboundary(fes2)
edge_fe_idx2 = selectelem(fens2, boundaryfes2, box=[0.5,0.5, 0.0,height2, 0.0,depth2], inflate=1e-8)
edge_fes2 = subset(boundaryfes2, edge_fe_idx2)
dbc_fes_idx2 = setdiff(1:count(boundaryfes2), edge_fe_idx2)

fens2.xyz[:, 1] .+= bend * (1.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 0.5).^2

geom2 = NodalField(fens2.xyz)
u2 = NodalField(zeros(size(fens2.xyz, 1), 3)) # displacement field

dbc_nodes2 = unique(stack(boundaryfes2.conn[dbc_fes_idx2]))
for i in dbc_nodes2
    setebc!(u2, [i], 1, ux(fens2.xyz[i, :]...))
    setebc!(u2, [i], 2, uy(fens2.xyz[i, :]...))
    setebc!(u2, [i], 3, uz(fens2.xyz[i, :]...))
end

applyebc!(u2)
numberdofs!(u2)
femm2 = FEMMDeforLinear(MR, IntegDomain(fes2, Rule2), material)
K2 = stiffness(femm2, geom2, u2)
K2_ff = matrix_blocked(K2, nfreedofs(u2), nfreedofs(u2))[:ff]
K2_fd = matrix_blocked(K2, nfreedofs(u2), nfreedofs(u2))[:fd]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(u2))[:f] - K2_fd * gathersysvec(u2, :d)

edge_nodes2 = selectnode(fens2; box=[1.0, 1.0, 0.0, height2, 0.0, depth2], inflate=1e-8)
##########################################################################################

xs_i = 0.5
ys_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
zs_i = collect(linearspace(0.0, 1.0, N_elem_i+1))
# fens_i, fes_i = T3blockx(ys_i, zs_i, :a)
fens_i, fes_i = Q4blockx(ys_i, zs_i)
fens_i.xyz = hcat(xs_i*ones(size(fens_i.xyz, 1), 1), fens_i.xyz)
fens_i.xyz[:, 1] .+= bend * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 0.5).^2

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 3)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 3  )) # Lagrange multipliers field
end
numberdofs!(u_i)
femm_i = FEMMDeforLinear(MR, IntegDomain(fes_i, Rule1), material)
@time D1, meta1 = common_refinement(fens1, edge_fes1, fens_i, fes_i; lam_order=lam_order, h=0.03, dim_u=3, tri_order = 2,)
@time D2, meta2 = common_refinement(fens2, edge_fes2, fens_i, fes_i; lam_order=lam_order, h=0.3, dim_u=3, tri_order = 2,)
# error("The 3D patch test is not implemented yet. Please use the 2D version instead.")

# D1,Pi_NC1,Pi_phi1 = build_D_matrix(fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
# D2,Pi_NC2,Pi_phi2 = build_D_matrix(fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)

dbc_dofs1 = sort([3*dbc_nodes1.-2; 3*dbc_nodes1.-1; 3*dbc_nodes1])
dbc_dofs2 = sort([3*dbc_nodes2.-2; 3*dbc_nodes2.-1; 3*dbc_nodes2])

dbc_lam_f = -(D1[:,dbc_dofs1] * gathersysvec(u1, :d) - D2[:,dbc_dofs2] * gathersysvec(u2, :d))
D1 = D1[:,setdiff(1:3*count(fens1), dbc_dofs1)]
D2 = D2[:,setdiff(1:3*count(fens2), dbc_dofs2)]

A = [K1_ff           spzeros(size(K1_ff,1), size(K2_ff,2))    D1';
     spzeros(size(K2_ff,1), size(K1_ff,2))    K2_ff           -D2';
     D1                     -D2                 spzeros(size(D1,1), size(D1,1))]
B = [F1_ff;
     F2_ff;
     dbc_lam_f]
X = A\B

scattersysvec!(u1, X[1:size(K1_ff,1)])
scattersysvec!(u2, X[size(K1_ff,1)+1:size(K1_ff,1)+size(K2_ff,1)])
scattersysvec!(u_i, X[size(K1_ff,1)+size(K2_ff,1)+1:end])

err1 = L2_err3D(femm, geom1, u1, ux, uy, uz)
err2 = L2_err3D(femm2, geom2, u2, ux, uy, uz)


# File1 = "Patch_1.vtk"
# vtkexportmesh(
#     File1,
#     fes1.conn,
#     geom1.values,
#     exporter1;
#     vectors = [("u", u1.values)],
# )
# File2 = "Patch_2.vtk"
# vtkexportmesh(
#     File2,
#     fes2.conn,
#     geom2.values,
#     exporter2;
#     vectors = [("u", u2.values)],
# )

filename = "3D_patch1.vtk"
vtkexportmesh(
    filename,
    fens1,
    fes1;
    scalars = [
        
        ("Err", err1.values)
    ],
    vectors = [
        ("u", u1.values)
    ]
)

filename = "3D_patch2.vtk"
vtkexportmesh(
    filename,
    fens2,
    fes2;
    scalars = [
        ("Err", err2.values)
    ],
    vectors = [
        ("u", u2.values)
    ]
)
# println(u_i.values)