using SparseArrays

function build_D_matrix(fens_i, fens_sd, boundarys_fes_sd; tol=1e-8)
    edge_nodes_sd = unique(vcat(boundarys_fes_sd...))
    fens_u, fes_u, M_u = build_union_mesh(fens_i, fens_sd, edge_nodes_sd)
    X = fens_u.xyz[ :, 1:2]
    Pi_phi = Lagrange_interpolation_matrix(X, fens_i.xyz[:, 1:2], fes_i.conn, 1; tol=tol)
    Pi_NC = Lagrange_interpolation_matrix(X, fens_sd.xyz[edge_nodes_sd, 1:2], boundarys_fes_sd.conn, 1; tol=tol)
    C = Pi_Phi.T * M_u * Pi_NC
    
end

function build_union_mesh(fens_i, fens_sd, edge_nodes_sd)
    endpoints = unique(vcat(fens_i.xyz[:, 1], fens_sd.xyz[edge_nodes_sd, 1]))
    fens_u, fes_u = L2blockx2D(endpoints, fens_i.xyz[:, 2])
    E = 1.0
    nu = 0.0
    rho = 1.0
    MR = DeforModelRed1D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    geom_u = NodalField(fens_u.xyz)
    u_u = NodalField(zeros(size(fens_u.xyz, 1), 2))
    numberdofs!(u_u)
    femm_u = FEMMDeforLinear(MR, IntegDomain(fes_u, GaussRule(1, 2)), material)
    M_u = mass(femm_u, geom_u, u_u)
    return fens_u, fes_u, M_u
end

function Lagrange_interpolation_matrix(X, Y, conn, p; tol = 1e-8)
    eps_end = tol
    npts = size(X, 1)
    nnds = size(Y, 1)
    nels = size(conn, 1)
    I = Int[]
    J = Int[]
    V = Float64[]
    for r in 1:npts
        for elem_idx in 1:nels
            nodes = conn[elem_idx, :]
            in_elem, xi, dist = point_in_element((Float64(X[r,1]), Float64(X[r,2])), Y, nodes, p; tol=tol)
            in_elem||continue
            if abs(xi + 1.0) <= eps_end
                push!(I,r)
                push!(J,nodes[1])
                push!(V,1.0)
                break
            elseif abs(xi - 1.0) <= eps_end
                push!(I,r)
                push!(J,nodes[end])
                push!(V,1.0)
                break
            else
                N = lagrange_1d(xi, p)
                for a in 1:(p+1)
                    push!(I,r)
                    push!(J,nodes[a])
                    push!(V,N[a])
                end
                break
            end
        end
    end
    # if max value is >1, give error
    if maximum(V) > 1.0 + tol
        error("Lagrange interpolation matrix has values greater than 1.")
    end

    sparse(I, J, V, npts, nnds)
end

function lagrange_1d(xi, p)
    xi_n = range(-1.0, 1.0; length = p+1)
    N = ones(p+1)
    for a in 1:(p+1)
        for b in 1:(p+1)
            if a != b
                N[a] *= (xi - xi_n[b]) / (xi_n[a] - xi_n[b])
            end
        end
    end
    return N
end

function lagrange_shapes_deriv_1d(p::Int, xi::Float64)
    if p == 1
        # nodes at [-1, 1]: dN1/dxi = -1/2, dN2/dxi = 1/2
        return [-0.5, 0.5]
    elseif p == 2
        # nodes at [-1, 0, 1]:
        # N0 = xi*(xi-1)/2, N1 = 1 - xi^2, N2 = xi*(xi+1)/2
        # dN = [xi - 1/2, -2*xi, xi + 1/2]
        return [xi - 0.5, -2.0*xi, xi + 0.5]
    else
        xin = range(-1.0, 1.0, length=p+1)
        dN = zeros(p+1)
        @inbounds for a in 1:p+1
            xia = xin[a]
            s = 0.0
            for k in 1:p+1
                k == a && continue
                num = 1.0; den = 1.0
                for b in 1:p+1
                    b == a && continue
                    if b == k
                        den *= (xia - xin[b])
                    else
                        num *= (xi - xin[b])
                        den *= (xia - xin[b])
                    end
                end
                s += num/den
            end
            dN[a] = s
        end
        return dN
    end
end


function curve_map(xi, Y, nodes, p)
    N   = lagrange_1d(xi, p)
    dN  = lagrange_shapes_deriv_1d(p, xi)
    ddN = lagrange_shapes_dd_1d(p, xi)
    x = 0.0; y = 0.0
    dx = 0.0; dy = 0.0
    for a in 1:p+1
        j = nodes[a]
        Xja, Yja = Y[j,1], Y[j,2]
        Na, dNa = N[a], dN[a]
        x   += Na  * Xja;    y   += Na  * Yja
        dx  += dNa * Xja;    dy  += dNa * Yja
    end
    return x, y, dx, dy
end

function get_xi(P, Y, nodes, p;
                     xi0=0.0, tol=1e-12, maxiter=20,
                     alpha=1e-12, maxstep=0.5)
    xi = xi0
    for _ in 1:maxiter
        x, y, dx, dy = curve_map(xi, Y, nodes, p)
        rx, ry = x - P[1], y - P[2]
        dist = hypot(rx, ry)
        b = dx*rx + dy*ry
        if abs(b) <= tol
            return xi, dist, true
        end
        A = dx*dx + dy*dy
        dxi = - b / (A + alpha)
        if abs(dxi) > maxstep
            dxi = sign(dxi) * maxstep
        end
        xi = clamp(xi + dxi, -1.0, 1.0)
    end
    x, y, _, _= curve_map(xi, Y, nodes, p)
    dist = hypot(x - P[1], y - P[2])
    return xi, dist, false
end



function point_in_element(P, Y,
                          nodes, p; tol)
    xi, d, ok = get_xi(P, Y, nodes, p)
    return (ok && d <= tol && xi >= -1.0 - 1e-12 && xi <= 1.0 + 1e-12), xi, d
end

x = [0 0; 1 1; 2 4; 3 9; 4 16]
y = [0  0; 1 1; 2 4; 3 9; 4 16]
conn = [1 2; 2 3; 3 4; 4 5]
P = Lagrange_interpolation_matrix(x, y, conn, 1)