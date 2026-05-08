using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
include("utilities.jl")


println("Q4. Plane stress.")

E = 1.0
nu = 1/3
MR = DeforModelRed2DStress
material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

N_elem1 = 20
N_elem2 = 30
N_elem_i = min(N_elem1, N_elem2)
left_m = "q"
right_m = "t"
skew = 0.0
bend = 0.2
lam_order = 0


alpha, beta, gamma, delta, eta, phi =
    1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
ux(x, y) = alpha + beta * x + gamma * y
uy(x, y) = delta + eta * x + phi * y

#########################################################################################
width1 = 1.0
height1 = 2.0
if left_m == "t"
    fens1, fes1 = T6block(width1, height1, floor(Int, N_elem1/2), N_elem1)
    Rule1 = TriRule(9)
    exporter1 = FinEtools.MeshExportModule.VTK.T6
else
    xs1 = collect(linearspace(0.0, width1, floor(Int, N_elem1/2)+1))
    ys1 = collect(linearspace(0.0, height1, N_elem1+1))
    fens1, fes1 = Q9blockx(xs1, ys1)
    Rule1 = GaussRule(2,4)
    exporter1 = FinEtools.MeshExportModule.VTK.Q9
end

boundaryfes1 = meshboundary(fes1)
edge_fe_idx1 = selectelem(fens1, boundaryfes1,  box=[width1,width1, 0.0,height1], inflate=1e-8)
edge_fes1 = subset(boundaryfes1, edge_fe_idx1)
dbc_fes_idx1 = setdiff(1:count(boundaryfes1), edge_fe_idx1)
p = maximum(length.(edge_fes1.conn)) - 1

#########################################################################################
width2 = 1.0
height2 = 2.0
if right_m == "t"
    fens2, fes2 = T6block(width2, height2, floor(Int, N_elem2/2), N_elem2)
    Rule2 = TriRule(3)
    # shift the second mesh to the right by 1.0
    fens2.xyz[:, 1] .+= 1.0
    exporter2 = FinEtools.MeshExportModule.VTK.T6
else
    xs2 = collect(linearspace(1.0 ,1.0+ width2, floor(Int, N_elem2/2)+1))
    ys2 = collect(linearspace(0.0, height2, N_elem2+1))
    fens2, fes2 = Q9blockx(xs2, ys2)
    Rule2 = GaussRule(2,4)
    exporter2 = FinEtools.MeshExportModule.VTK.Q9
end

boundaryfes2 = meshboundary(fes2)
edge_fe_idx2 = selectelem(fens2, boundaryfes2, box=[1.0,1.0, 0.0,height2], inflate=1e-8)
edge_fes2 = subset(boundaryfes2, edge_fe_idx2)
dbc_fes_idx2 = setdiff(1:count(boundaryfes2), edge_fe_idx2)
#########################################################################################
xs_i = ones(N_elem_i+1)
ys_i = collect(linearspace(0.0, 2.0, N_elem_i+1))
fens_i, fes_i = L3blockx2D(xs_i, ys_i)

fens_u1, fes_u1, _ = build_union_mesh(fens_i,fes_i, fens1, edge_fes1, p; lam_order=lam_order)
fens_u2, fes_u2, _ = build_union_mesh(fens_i,fes_i, fens2, edge_fes2, p; lam_order=lam_order)
#############################################################################################


# fens1.xyz[:, 1] .+= skew * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 1.0)
fens1.xyz[:, 1] .+= bend * fens1.xyz[:, 1].*(fens1.xyz[:, 2] .- 1.0).^2


geom1 = NodalField(fens1.xyz)
u1 = NodalField(zeros(size(fens1.xyz, 1), 2)) # displacement field

dbc_nodes1 = unique(stack(boundaryfes1.conn[dbc_fes_idx1]))
for i in dbc_nodes1
    setebc!(u1, [i], 1, ux(fens1.xyz[i, :]...))
    setebc!(u1, [i], 2, uy(fens1.xyz[i, :]...))
end

applyebc!(u1)
numberdofs!(u1)
femm = FEMMDeforLinear(MR, IntegDomain(fes1, Rule1), material)

K1 = stiffness(femm, geom1, u1)
K1_ff = matrix_blocked(K1, nfreedofs(u1), nfreedofs(u1))[:ff]
K1_fd = matrix_blocked(K1, nfreedofs(u1), nfreedofs(u1))[:fd]
F1 = zeros(size(K1, 1))
F1_ff = vector_blocked(F1, nfreedofs(u1))[:f] - K1_fd * gathersysvec(u1, :d)

edge_nodes1 = selectnode(fens1; box=[width1,width1, 0.0,height1], inflate=1e-8)
#########################################################################################



fens2.xyz[:, 1] .+= bend * (2.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 1.0).^2

# fens2.xyz[:, 1] .+= skew * (2.0 .-fens2.xyz[:, 1]).*(fens2.xyz[:, 2] .- 1.0)

geom2 = NodalField(fens2.xyz)
u2 = NodalField(zeros(size(fens2.xyz, 1), 2)) # displacement field

dbc_nodes2 = unique(stack(boundaryfes2.conn[dbc_fes_idx2]))
for i in dbc_nodes2
    setebc!(u2, [i], 1, ux(fens2.xyz[i, :]...))
    setebc!(u2, [i], 2, uy(fens2.xyz[i, :]...))
end

applyebc!(u2)
numberdofs!(u2)
femm2 = FEMMDeforLinear(MR, IntegDomain(fes2, Rule2), material)
K2 = stiffness(femm2, geom2, u2)
K2_ff = matrix_blocked(K2, nfreedofs(u2), nfreedofs(u2))[:ff]
K2_fd = matrix_blocked(K2, nfreedofs(u2), nfreedofs(u2))[:fd]
F2 = zeros(size(K2, 1))
F2_ff = vector_blocked(F2, nfreedofs(u2))[:f] - K2_fd * gathersysvec(u2, :d)

edge_nodes2 = selectnode(fens2; box=[1.0, 1.0, 0.0, height2], inflate=1e-8)
##########################################################################################

# fens_i.xyz[:, 1] .+= skew * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 1.0)
fens_i.xyz[:, 1] .+= bend * fens_i.xyz[:, 1].*(fens_i.xyz[:, 2] .- 1.0).^2
fens_u1.xyz[:, 1] .+= bend * fens_u1.xyz[:, 1].*(fens_u1.xyz[:, 2] .- 1.0).^2
fens_u2.xyz[:, 1] .+= bend * fens_u2.xyz[:, 1].*(fens_u2.xyz[:, 2] .- 1.0).^2

geom_i = NodalField(fens_i.xyz)
if lam_order == 0
    u_i  = ElementalField(zeros(count(fes_i), 2)) # Lagrange multipliers field
else
    u_i  = NodalField(zeros(size(fens_i.xyz, 1), 2)) # Lagrange multipliers field
end
numberdofs!(u_i)
D1,Pi_NC1,Pi_phi1 = build_D_matrix(fens_u1, fes_u1, fens_i, fes_i, fens1, edge_fes1; lam_order=lam_order,tol=1e-8)
D2,Pi_NC2,Pi_phi2 = build_D_matrix(fens_u2, fes_u2, fens_i, fes_i, fens2, edge_fes2; lam_order=lam_order,tol=1e-8)

dbc_dofs1 = sort([2*dbc_nodes1.-1; 2*dbc_nodes1])
dbc_dofs2 = sort([2*dbc_nodes2.-1; 2*dbc_nodes2])

dbc_lam_f = -(D1[:,dbc_dofs1] * gathersysvec(u1, :d) - D2[:,dbc_dofs2] * gathersysvec(u2, :d))
D1 = D1[:,setdiff(1:2*count(fens1), dbc_dofs1)]
D2 = D2[:,setdiff(1:2*count(fens2), dbc_dofs2)]

A = [K1_ff           zeros(size(K1_ff,1), size(K2_ff,2))    D1';
     zeros(size(K2_ff,1), size(K1_ff,2))    K2_ff           -D2';
     D1                     -D2                 zeros(size(D1,1), size(D1,1))]
B = [F1_ff;
     F2_ff;
     dbc_lam_f]
X = A\B

scattersysvec!(u1, X[1:size(K1_ff,1)])
scattersysvec!(u2, X[size(K1_ff,1)+1:size(K1_ff,1)+size(K2_ff,1)])
scattersysvec!(u_i, X[size(K1_ff,1)+size(K2_ff,1)+1:end])


err1 = L2_err(femm, geom1, u1, ux, uy)
err2 = L2_err(femm2, geom2, u2, ux, uy)

filename = "curved1.vtk"
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

filename = "curved2.vtk"
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