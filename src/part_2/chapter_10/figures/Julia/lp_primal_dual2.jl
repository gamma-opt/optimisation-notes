################################################################################
## Test code for primal-dual IP                                               ##
################################################################################
using JuMP
using Clp
using Plots
using LinearAlgebra
using LaTeXStrings

pyplot()


function strict_primal(c::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})

    (m,n) = size(A)

    lp = Model(with_optimizer(Clp.Optimizer))

    @variable(lp, t >= 0)
    @variable(lp, x[1:n] >= 0)
    @constraint(lp, A*x .== b)
    @constraint(lp, [i = 1:n], t >= 1 - x[i])
    @objective(lp, Min, t)

    optimize!(lp)

    return (value(t), value.(x))
end

## Update Newton direction
function update_newton_dir(c::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64},
                           x::Vector{Float64}, u::Vector{Float64}, v::Vector{Float64},
                           μ::Float64)
    ## Diagonalize u and x
    U = Diagonal(u)
    X = Diagonal(x)

    ## Calculate parmeters
    e = ones(length(x))
    (m,n) = size(A)

    ## Update Newton direction
    Nmatrix = [A zeros(m,m) zeros(m,n);
           zeros(n,n) A' Diagonal(ones(n));
           U zeros(n,m) X]

    rp = A*x .- b
    rd = A'v .+ u .- c
    rc = X*U*e - μ*e

    Nrhs = -[rp; rd; rc]

    d = Nmatrix \ Nrhs

    dₓ = d[1:n]
    dᵥ = d[n+1:n+m]
    dᵤ = d[n+m+1:end]
    return (dᵥ,dᵤ,dₓ)
end

################################################################################
## Calculate step size that retains feasibility                             ##
################################################################################
function calculate_step_size(x::Vector{Float64}, d::Vector{Float64}, ϵ::Float64)
    n = length(d)
    α = 0.9999
    for i=1:n
        if d[i] < 0
            α = min(α, -x[i]/d[i]) # prevents variable becoming negative
        end
    end
    return round(α, digits=Int(-log10(ϵ))) #rounding avoids numerical issues
end

## For validating the optimal solution
function linprog(c, A, b)

    m,n = size(A)

    lp = Model(with_optimizer(Clp.Optimizer))

    @variable(lp, x[1:n])
    @constraint(lp, cons, A*x .== b)
    @constraint(lp, bound, x .>= 0)
    @objective(lp, Min, sum(c[i]*x[i] for i=1:n))

    optimize!(lp)

    return (value.(x), dual.(bound), dual.(cons))
end

function primal_dual_ip(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64},
                        β::Float64, ϵ::Float64, N::Int)

################################################################################
## Solve the LP problem and its dual                                          ##
##                                                                            ##
## minimize    cᵀx          maximize    bᵀv                                   ##
## subject to  Ax = b       subject to  Aᵀv + u = c                           ##
##              x ≥ 0                         u ≥ 0                           ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
################################################################################

    n = length(c)    # Number of primal variables x and dual variables u
    m = length(b)    # Number of dual variables v
    x = zeros(N, n)  # Primal variable x values
    u = zeros(N, n)  # Dual variable u values
    v = zeros(N, m)  # Dual variable v values
    e = ones(n)      # Vector of ones
    α_p = 1          # Step size for variable x update
    α_d = 1          # Step size for variable u, and v update
    μ = 100.0         # Around 10 for efficiency.
                     # Large values highlight central path

    ## Find an initial solution
    (⋅, x₀) = strict_primal(c, A, b)
    u₀ = μ ./ x₀

    u[1,:] = u₀
    x[1,:] = x₀ # Note v₀ initialised as 0

    ## Main loop
    for i = 1:N
        ## Stopping condition #1
        if n*μ < ϵ           ## if dot(c,x[i,:]) - dot(b,v[i,:]) < ϵ
            v = v[1:i,:]
            u = u[1:i,:]
            x = x[1:i,:]
            return (v,u,x)
        end
        ## Update dₓ, dᵥ, dᵤ
        (dᵥ,dᵤ,dₓ) = update_newton_dir(c,A,b,x[i,:],u[i,:],v[i,:],μ)

        ## Calculate step size α primal and α dual
        α_p = calculate_step_size(x[i,:], dₓ, ϵ)
        α_d = calculate_step_size(u[i,:], dᵤ, ϵ)

        ## Update variables
        v[i+1,:] = v[i,:] + α_d*dᵥ
        u[i+1,:] = u[i,:] + α_d*dᵤ
        x[i+1,:] = x[i,:] + α_p*dₓ

        ## Stopping condition #2
        if i == N-1
            return (v,u,x)
        end
        ## Update μ
        μ = μ*β
    end
end



## Problem in canonical form (notice that it's a 2D problem)
# min      -x₁ - x₂
# s.t.: -1/3x₁ + x₂ ≦  5
#        1/5x₁ - x₂ ≦ -1
#       -8/3x₁ - x₂ ≦ -8
#        1/2x₁ + x₂ ≦  9
#           x₁ - x₂ ≦  4
#           x₁  ,x₂ ≧  0

## Problem data in standard form
c = [-1.0,-1.0, 0.0,0.0,0.0,0.0,0.0]
b = [ 5.0,-1.0,-8.0,9.0,4.0]
A = [-1/3  1.0 1.0 0.0 0.0 0.0 0.0;
      1/5 -1.0 0.0 1.0 0.0 0.0 0.0;
     -8/3 -1.0 0.0 0.0 1.0 0.0 0.0;
      1/2  1.0 0.0 0.0 0.0 1.0 0.0;
      1.0 -1.0 0.0 0.0 0.0 0.0 1.0]


# c = [2.0, 2.0, 0.0, 0.0]
# b = [-8.0; -10.0 ]
# A = [ -2.0  -1.0 1.0 0.0 ;
#     -1.0  -2.0 0.0 1.0 ]


# x = [10, 10]
# x₀ = [x[1], x[2], b[1] - A[1,1:2]'x, b[2] - A[2,1:2]'x]


## Executions and plotting
N = 500     # Number of iterations
β = 0.3     # Reduction factor
ϵ = 1e-6    # Tolerance

## Comparing solution time and the solution itself
@time (v,u,x) = primal_dual_ip(A, b, c, β, ϵ, N)
@time (x_opt, u_opt, v_opt) = linprog(c, A, b)

norm(x[end,:] - x_opt) ## This should be close to 0

## Plotting...
x1 = range(0,stop=15,length=100)
x2 = range(0,stop=15, length=100)

# plot(x1, 8 .- 2*x1,  color = :1, reuse = false)
# plot!(x1, (10 .- x1)./2,      color = :1,
# xaxis = ("x₁", (0,15)),
# yaxis = ("x₂", (0,15)),
# aspect_ratio = :equal,
# #size = (1000,800),
# legend = false,
# #title = "Optimization step"
# )

# annotate!([(0.2, 8, text(L"$2x_1^2 + x_2 \geq 8$",9,:left,
#                RGBA(0.0,0.605603,0.97868,1.0)))])
#
#
# annotate!([(8, 1.3, text(L"$x_1^2 + 2x_2 \geq 10$",9,:left,
#                RGBA(0.0,0.605603,0.97868,1.0)))])

plot(x1, 5 .+ (1/3)*x1,  color = :1, reuse = false)
plot!(x1, 1 .+ (1/5)*x1, color = :1)
plot!(x1, 8 .- (8/3)*x1, color = :1)
plot!(x1, 9 .- (1/2)*x1, color = :1)
plot!(x1, -4 .+ x1,      color = :1,
    xaxis = ("x₁", (0,10)),
    yaxis = ("x₂", (0,10)),
    aspect_ratio = :equal,
    #size = (800,800),
    legend = false,
    title = L"$\beta = 0.3, \mu = 100$"
    )

traj = x[1:end,1:2]
plot!(traj[:,1], traj[:,2], marker = :o, color = :3)
annotate!([(traj[1,1], traj[1,2]+0.1, text(L"$x_0$",9,:bottom,
                RGBA(0.2422242978521988,0.6432750931576304,0.30444865153411527,1.0) ))])

f(x) = x[1] + x[2]
f(traj[end,:])

plot!(x1, f(traj[end,:]) .- x1, color = :2, line = :dash)

# annotate!([(10, 9.8, text(L"$x_0$",9,:top,
#                 RGBA(0.2422242978521988,0.6432750931576304,0.30444865153411527,1.0) ))])

savefig("lp_example2.pdf")

# Plot function contours
# contour!(x,y,z, # Change title accordingly
#         levels = [0.25, 1, 2.5 , 5, 9],
#         xaxis = (L"$x_1$", (0,2)),
#         yaxis = (L"$x_2$", (0,2)),
#         clims = (0,12),
#         clabels = true,
#         legend = false,
#         aspect_ratio = :equal)
