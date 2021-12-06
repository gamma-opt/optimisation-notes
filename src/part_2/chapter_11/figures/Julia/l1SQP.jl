using JuMP, Ipopt, ForwardDiff, LinearAlgebra

∇(f,x) = ForwardDiff.gradient(f,x)
∇²(f,x) = ForwardDiff.hessian(f,x)

function l1SQP()
        # Problem information
        n = 2 # total of variables
        m = 4 # total of ineq. constraints

        f(x) =  2x[1]^2 + 2x[2]^2 - 2(x[1]*x[2]) - 4x[1] - 6x[2]
        ∇f(x) = ∇(f,x)
        ∇²f(x) = ∇²(f,x)

        g₁(x) = x[1]^2 - x[2]   # ≦ 0
        g₂(x) = x[1] + 5x[2] - 5 # ≦ 0
        g₃(x) = -x[1] # ≦ 0
        g₄(x) = -x[2] # ≦ 0

        ∇g₁(x) = ∇(g₁,x)
        ∇²g₁(x) = ∇²(g₁,x)
        ∇g₂(x) = ∇(g₂,x)
        ∇g₃(x) = ∇(g₃,x)
        ∇g₄(x) = ∇(g₄,x)

        L(x,u) = f(x) + u[1]*g₁(x) + u[2]*g₂(x) + u[3]*g₃(x) + u[4]*g₄(x)
        ∇²L(x,u) = ∇²f(x) + u[1]*∇²g₁(x)  # g₂, g₃, and g₄ vanish as they are linear .

# Initialisation

        k = 1             # iter count
        N = 10            # max iterations
        ϵ = 1e-6          # tolerance
        μ = 10            # penalty
        Δᵏ = 5            # trust region
        xᵏ = [0.5; 0.5]     # starting point
        uᵏ = zeros(4)     # initial dual values
        x = zeros(2, N)   # trajectory
        x[:,k] = xᵏ

        ## Main loop
        for k = 1:N-1
                xᵏ = x[:,k]
                # Calculate ∇f, ∇h, and ∇L²
                ∇fᵏ = ∇f(xᵏ)
                ∇gᵏ = [∇g₁(xᵏ), ∇g₂(xᵏ), ∇g₃(xᵏ), ∇g₄(xᵏ)]
                ∇L²ᵏ = ∇²L(xᵏ,uᵏ)
                gᵏ = [g₁(xᵏ), g₂(xᵏ), g₃(xᵏ), g₄(xᵏ)]

                # Direction search subproblem
                QP = Model()
                set_optimizer(QP, Ipopt.Optimizer)
                @variable(QP, d[1:n])
                @objective(QP, Min, dot(∇fᵏ,d) + 0.5*(d'*∇L²ᵏ*d))
                @constraint(QP, LinearIneq[i=1:m], dot(∇gᵏ[i],d) + gᵏ[i] <= 0)

                optimize!(QP)

                dᵏ = value.(d)           # obtain direction
                x[:,k+1] = xᵏ + dᵏ       # step
                uᵏ = dual.(LinearIneq)   # obtain dual values

                if norm(x[:,k+1] - x[:,k]) < ϵ # stop condition
                        x = x[:,1:k+1]
                        return x
                end
        end
end

f(x) =  2x[1]^2 + 2x[2]^2 - 2(x[1]*x[2]) - 4x[1] - 6x[2]
x = l1SQP()

## Plotting...
n = 1000
x1 = range(0,5,length=n)
x2 = range(0,5,length=n)
z = [f([x1[i],x2[j]]) for j=1:n, i=1:n]

using Plots, LaTeXStrings
pyplot()

contour(x1, x2, z,
        clims = (-10,-5))
        #clabels = true)
plot!(x1, x1.^2, fill=(10,0.1), color=3)
plot!(x1, 1 .- (1/5)*x1, fill=(0,0.1), color=2,
    xaxis = (L"$x_1$", (0,2)),
    yaxis = (L"$x_2$", (0,2)),
    aspect_ratio = :equal,
    colorbar = false,
    legend = false)

plot!(x[1,:], x[2,:], marker = :o, color=1)

# First order approximation at initial point for constraint 1
plot!(x1, -0.25 .+ x1, color=3, line = :dash) # dot(∇gᵏ[i],d) + gᵏ[i] = 0


annotate!([(0.5, 0.4, text(L"$x_0$",9,:bottom,
                RGBA(0.0,0.605603,0.97868,1.0) ))])

annotate!([(0.55, 0.25, text(L"$g(x_0) + \nabla g(x_0)^\top x$",9,:left,
                color = :green))])


#savefig("l1sqp.pdf")
