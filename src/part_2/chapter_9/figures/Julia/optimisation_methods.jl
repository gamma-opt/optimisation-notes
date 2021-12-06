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
