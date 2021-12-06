using JuMP
using Cbc, Ipopt
using Random
using LinearAlgebra
using ForwardDiff
using Plots, LaTeXStrings
pyplot()

∇(f, x) = ForwardDiff.gradient(f, x)
D(θ, λ) = ForwardDiff.derivative(θ, λ)

function Armijo_ls(θ, λ, α, β, λbar)
    θ₀  = θ(0)                                 # Function value at zero
    Dθ₀ = D(θ, 0)                              # Derivative (slope) at zero
    while (θ(λ) > θ₀ + α*λ*Dθ₀) && (λ <= λbar) # Check termination condition
        λ = β*λ                                # Reduce λ until condition is satisfied
    end
    return λ
end

function FW()

    ## Constants
    λbar = 1.0            # Upper bound of λ in FW method
    λ    = λbar           # Initial step size λ
    eps  = 1e-4          # Convergence tolerance

    ## Armijo parameters
    α   = 0.01
    β   = 0.70

    ## Objective function to be minimized
    f(x) = exp(-(x[1]-3)/2) + exp((4x[2] + x[1] - 20)/10) + exp((-4x[2] + x[1])/10)
    ∇f(x) = ∇(f,x)
    A = [2 3;
         1 4]
    b = [8; 6]

    ## Initial values
    xᵏ  = [0; 0]
    traj = zeros(2,50)
    ## Compute x̄ᵏ from x̄ᵏ = argmin{x ∈ ℜⁿ : ∇f(xᵏ)ᵀx, x ∈ S}
    model = Model()
    set_optimizer(model, Cbc.Optimizer)
    @variable(model, x[1:2] >= 0)
    @constraint(model, A * x .<= b)
    @objective(model, Min, dot(∇(f, xᵏ), (x - xᵏ)))
    optimize!(model)
    x̄ᵏ = value.(x)
    dᵏ = x̄ᵏ - xᵏ

    ## Iteration counter + 1st solution
    k   = 1
    xᵏ  = xᵏ + λ*dᵏ
    traj[:,k+1] = xᵏ
    while abs(dot(∇(f, xᵏ), dᵏ)) > eps

        ## Compute x̄ᵏ from x̄ᵏ = argmin{x ∈ ℜⁿ : ∇f(xᵏ)ᵀx, x ∈ S}
        model = Model()
        set_optimizer(model, Cbc.Optimizer)
        @variable(model, x[1:2] >= 0)
        @constraint(model, A * x .<= b)
        @objective(model, Min, dot(∇(f, xᵏ), (x - xᵏ)))
        optimize!(model)
        x̄ᵏ = value.(x)
        traj[:,k+1] = xᵏ
        dᵏ = x̄ᵏ - xᵏ

        #### Line search
        θ(λ) = f(xᵏ + λ*dᵏ)
        λ    = Armijo_ls(θ, λ, α, β, λbar)

        #### Update solution
        k  = k + 1
        xᵏ = xᵏ + λ*dᵏ
        println("residual: ", round(dot(∇(f, xᵏ), dᵏ), digits=5), " / iter: ", k)
    end

    ## Get optimal cost
    obj = f(xᵏ)
    return (k, xᵏ, obj, traj[:,1:k])
end

## Solve model
(k, xsol, obj, traj) = FW()

## Print solution
println("  Iterations: ", k)
println("Optimal cost: ", round(obj, digits = 2))

# Plotting

# Plotting the contours of the function optimised and the trajectory

f(x) = exp(-(x[1]-3)/2) + exp((4x[2] + x[1] - 20)/10) + exp((-4x[2] + x[1])/10)
# For validation
A = [2 3;
     1 4]
b = [8; 6]


n = 1000
x = range(-20,stop=20,length=n);
y = range(-20,stop=20,length=n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

# Plot function contours
contour(x,y,z, # Change title accordingly
        levels = [0.5, 1, 1.5, 2, 2.5, 3, 4],
        xaxis = (L"$x_1$", (-0.1,6)),
        yaxis = (L"$x_2$", (-0.1,5)),
        clims = (0,5),
        clabels = true,
        legend = false)
        #aspect_ratio = :equal)

# Plot feasible region
plot!(x, ((b[1] .- A[1,1] * x)./A[1,2]),
        #fill= (10,0.1),
        color = 1)

plot!(x, ((b[2] .- A[2,1] * x)./A[2,2]),
        #fill= (10,0.1),
        color = 2)
#annotate!([(0.25, 1.75, text(L"$x_1^2 - x_2 \leq 0$",9,:left,
#                RGBA(0.0,0.605603,0.97868,1.0)))])

# Plot trajectory
plot!(traj[1,1:3], traj[2,1:3], marker=:o, color = 3)
annotate!([(0.1, 0.1, text(L"$x_0$",9,:left,
                 RGBA(0.2422242978521988,0.6432750931576304,0.30444865153411527,1.0) ))])

savefig("FW_first.pdf")

plot!(traj[1,:], traj[2,:], marker=:o, color = 3)

savefig("FW_alliter.pdf")



# For validation
A = [2 3;
     1 4]
b = [8; 6]

model = Model(with_optimizer(Ipopt.Optimizer))
@variable(model, x[1:2] >= 0)
@constraint(model, A * x .<= b)
@NLobjective(model, Min, exp(-(x[1]-3)/2) + exp((4x[2] + x[1] - 20)/10) + exp((-4x[2] + x[1])/10))
optimize!(model)
value.(x)
