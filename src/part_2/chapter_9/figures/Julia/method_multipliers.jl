using JuMP, Ipopt, LinearAlgebra
using ForwardDiff

# Line search
function bisection(f,a,b)
    l = 1e-12
    bᵢ = b
    aᵢ = a
    λ = 0
    while bᵢ - aᵢ > l
        λ = (aᵢ + bᵢ)/2
        ∇f= ForwardDiff.derivative(f,λ)
        if abs(∇f) < l # equivalent to ∇f = 0
            return λ
        elseif ∇f > 0  # left move
            aᵢ₊₁ = aᵢ
            bᵢ₊₁ = λ
        else           # move right
            aᵢ₊₁ = λ
            bᵢ₊₁ = bᵢ
        end
        aᵢ = aᵢ₊₁
        bᵢ = bᵢ₊₁
    end
        return λ
end;

function newton(f, x, N = 1000; ϵ = 1e-6, λ = 1.0, a = 0.0, b = 2.0)
    tini = time()                  # Start timer
    ∇(f, x)  = ForwardDiff.gradient(f, x)
    H(f, x) = ForwardDiff.hessian(f, x)

    for i = 1:N
        ∇fᵢ = ∇(f, x)               # Gradient
        if norm(∇fᵢ) < ϵ            # Stopping condition #1
            tend = time() - tini    # Computation time
            return x, f(x), tend, i # Return cost, time, iterations
        end
        d = -H(f, x)\∇fᵢ            # Newton direction
        ls(λ) = f(x + λ*d)          # Line search via bisection
        λ = bisection(ls,a,b)
        x = x + λ*d                 # Move to a new point
    end
    tend = time() - tini            # Computation time
    return x, f(x), tend, N         # Return x, f(x), time, iterations
end;

# We will consider this problem
f(x) = x[1]^2+ x[2]^2   # min f(x)
h(x) = x[1] + x[2] - 1  # h(x) = 0

# First, we will solve it using JuMP and Ipopt
m = Model()
@variable(m, x[1:2] >= 0)
@objective(m, Min, sum(x[i]^2 for i = 1:2))
@constraint(m, sum(x[i] for i = 1:2) == 1)

optimize!(m, with_optimizer(Ipopt.Optimizer))
println("Optimal solution from Ipopt: ", value.(x))

# Now, let's see if we can obtain the same results with
# an implementation of the augmented Lagrangian method of multipliers

# Initial parameters
μ = 10      # penalty term. Try different values
ϵ = 1e-6    # tolerance
xᵏ = [0,0]  # initial point
vᵏ = 0      # initial dual value

# Method of multipliers implementation
function method_of_multipliers(f, h, xᵏ = [0,0], vᵏ = 0, μ = 10, ϵ = 1e-6)
    println("\nStarting method of multipliers...")
    println("Starting point: ", xᵏ, "/ Objective : ", f(xᵏ))
    k = 0
    while (abs(h(xᵏ)) > ϵ)
        L(x) = f(x) + vᵏ*(h(x)) + μ*(h(x))^2 # aug. Lag. function
        xᵏ,⋅,⋅,⋅ = newton(L,xᵏ)              # primal step using Newton's
        println("Current point: ", xᵏ, "/ Objective : ", f(xᵏ))
        vᵏ = vᵏ + 2μ*(h(xᵏ))                 # dual step
        k += 1
    end
    println("Converged in ",k , " iterations. Optimal: ", xᵏ)
    return xᵏ
end

method_of_multipliers(f,h)

f(x) = x[1]^2+ x[2]^2   # min f(x)
h(x) = x[1] + x[2] - 1  # h(x) = 0

# Try a different problem executing these:
f(x) = (x[1]-2)^4+ (x[1] - 2x[2])^2   # min f(x)
h(x) = x[1]^2 - x[2]                  # h(x) = 0

method_of_multipliers(f,h,[2,1])
