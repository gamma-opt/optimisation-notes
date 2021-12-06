using LinearAlgebra

# This helper function allows us changing the line search method used
# Creating control parameter for line searches

include("line_searches.jl");

# Calculates first and second-order derivatives
∇(f,x) = ForwardDiff.gradient(f, x);
H(f,x) = ForwardDiff.hessian(f, x);
ϵ_standard = 1e-4

function line_search(f, method)
    λbar = 0.0
    max_step = 5.0
    if method == :armijo
        λbar = armijo(f)
    elseif method == :bisection
        λbar = bisection(f, 0.0, max_step) 
    elseif method == :golden_ratio 
        λbar = golden_ratio(f, 0.0, max_step)
    elseif method == :dichotomous
        λbar = dichotomous(f, 0.0, max_step)    
    elseif method == :fibonacci
        λbar = fibonacci(f, 0.0, max_step)     
    elseif method == :newton
        λbar = newton(f, 1.0) # λ₀=1 prevents divergence 
    elseif method == :uniform
        λbar = uniform(f, 0.0, max_step, 5, 10)        
    elseif method == :fixed
        λbar = 1.0
    else
        println("Invalid line search method")
    end

    return λbar
end

# Method 1: coordinate descent
function coordinate_descent(f, xstart, M=1001, ϵ=ϵ_standard; verbose=false)
    n = length(xstart)
    xc = zeros(n,M) # store step history
    xc[:,1] = xstart # initial point
    k = 1; # counter

    println("Starting coordinate descent...")
    tini = time(); # Start stopwatch
    while k <= M-1
        # stop criteria on the norm of the difference between 2 sucessive points
        if   k > 1 && norm(xc[:,k] - xc[:,k-1]) < ϵ
             xc = xc[:,1:k]
             break
        end

        for j = 1:n
            d = zeros(n) # setting search direction as dⱼ = 1.
            d[j] = 1

            ls(λ) = f(xc[:,k] + λ*d) # state line search (ls) function
            λbar = dichotomous(ls, -5, 5) # Notice it requires λ ∈ [-5, 5]
            xc[:,k+1] = xc[:,k] + λbar*d

            if verbose
                println("iter=", k, " λ=", round(λbar, digits=2), " xᵏ=", xc[:,k+1])
            end

            k = k+1
        end
    end
    tend = time() - tini; # Stop stopwatch

    # Reporting solution found
    println("Coordinate descent converged.")
    println(" Total steps: ", k-1)
    println(" Total time (s): ", tend)
    println(" Sol. found: ", xc[:,k], "/ Opt. value: ", f(xc[:,k]), "\n")

    return xc
end

# Method 2: gradient descent
function gradient(f, xstart, ls_method=:dichotomous; M=1001, ϵ=ϵ_standard)
    n = length(xstart)
    xg = zeros(n,M) # store step history
    xg[:,1] = xstart # initial point
    k = 1

    println("Starting gradient descent...")
    tini = time(); # Start stopwatch
    while k <= M-1
        ∇f = ∇(f,xg[:,k])
        d = -∇f/norm(∇f)
        # stop criteria on the norm of ∇f(x)
        if norm(∇f) < ϵ
             xg = xg[:,1:k]
             break
        end
        # state line search (ls) function
        ls(λ) = f(xg[:,k] + λ*d)

        # Line search according to selected method
        λbar = line_search(ls, ls_method)

        xg[:,k+1] = xg[:,k] + λbar*d

        # uncomment to see progress of the algorithm
        #println("iter=", k, " λ=",λbar, " xᵏ=",xg[:,k+1])
        k = k+1
    end

    tend = time() - tini; # Stop stopwatch
    println("Gradient descent converged.")
    println(" Total steps: ", size(xg,2)-1)
    println(" Total time (s): ", tend)
    println(" Sol. found: ", xg[:,k], "/ Opt. value: ", f(xg[:,k]), "\n")

    return xg
end

# Method 3: Newton's method
function newton_method(f, xstart, ls_method=:dichotomous; M=1001, ϵ=ϵ_standard)
    n = length(xstart)
    xn = zeros(n,M) # store step history
    xn[:,1] = xstart # initial point
    k = 1

    println("Starting Newton's method...")
    tini = time(); # Start stopwatch
    while k <= M-1
        d = -H(f,xn[:,k])\∇(f,xn[:,k])
        #stop criteria on the norm of ∇f(x)
        if norm(d) < ϵ
             xn= xn[:,1:k]
             break
        end
        #state line search (ls) function
        ls(λ) = f(xn[:,k] + λ*d)

        # Line search according to selected method
        λbar = line_search(ls, ls_method)

        xn[:,k+1] = xn[:,k] + λbar*d
        # Uncomment to see progress of the algorithm
        #println("iter.=", k, " λ=",λbar, "xᵏ=",xg[:,k+1])

        k = k+1
    end

    tend = time() - tini; # Stop stopwatch
    println("Newton's method converged.")
    println(" Total steps: ", size(xn,2)-1)
    println(" Total time (s): ", tend)
    println(" Sol. found: ", xn[:,k], "/ Opt. value: ", f(xn[:,k]), "\n")

    return xn
end


### Lecture 6 methods


#Method 4: Conjugate gradient method
function conjugate_gradient(f, xstart, ls_method=:dichotomous; M=1001, ϵ=ϵ_standard)
    n = length(xstart)
    xj = zeros(n,M)
    xj[:,1] = xstart
    α = 0
    d = -∇(f,xstart)
    k = 1

    println("Starting conjugate gradient method...")
    tini = time(); # Start stopwatch
    while k <= M-1
        for j = 1:length(xstart) #n=2
            #state line search (ls) function
            ls(λ) = f(xj[:,k] + λ*d)

            # Line search according to selected method
            #λbar = line_search(ls, ARMIJO)
            λbar = line_search(ls, ls_method)

            xj[:,k+1] = xj[:,k] + λbar*d
            # uncomment to see progress of the algorithm
            #println("iter=", k, " λ=",λbar, " xᵏ=",xj[:,k+1])

            #α-update following Fletcher-Reeves
            α = norm(∇(f,xj[:,k+1]))^2/norm(∇(f,xj[:,k]))^2
            d = -∇(f,xj[:,k+1]) + α*d
            k = k+1
        end
        d = -∇(f,xj[:,k])
        # stop criteria on the norm of ∇f(x)
        if norm(d) < ϵ
             xj = xj[:,1:k]
             break
        end
    end
    tend = time() - tini; # Stop stopwatch
    println("Conjugate gradient converged.")
    println(" Total steps: ", size(xj,2)-1)
    println(" Total time (s): ", tend)
    println(" Sol. found: ", xj[:,k], "/ Opt. value: ", f(xj[:,k]), "\n")

    return xj
end

#Method 5: Quasi-Newton method (BFGS)
function bfgs(f, xstart, ls_method=:dichotomous; M=1001, ϵ=ϵ_standard)
    n = length(xstart)
    xb = zeros(n,M)
    xb[:,1] = xstart
    B = I(n)
    k = 1

    println("Starting quasi-Newton (BFGS) method...")
    tini = time(); # Start stopwatch
    while k <= M-1
        d = - B \ ∇(f, xb[:,k])

        ls(λ) = f(xb[:,k] + λ*d)

        # Line search according to selected method
        λbar = line_search(ls, ls_method)

        # Update function difference and do a step
        p = λbar*d
        xb[:,k+1] = xb[:,k] + p

        # Check stop criterion using the norm of ∇f(x)
        if norm(∇(f, xb[:,k+1])) < ϵ
            xb = xb[:,1:k+1]
            break
        end

        # Update Gradient difference
        q = ∇(f, xb[:,k+1]) - ∇(f, xb[:,k])
        # Update Hessian approximation
        B = B + (q*q')/(q'*p) - (B*p*p'*B)/(p'*B*p)

        k = k+1
    end
    tend = time() - tini; # Stop stopwatch
    println("Quasi-Newton (BFGS) converged.")
    println(" Total steps: ", size(xb,2)-1)
    println(" Total time (s): ", tend)
    println(" Sol. found: ", xb[:,k], "/ Opt. value: ", f(xb[:,k]), "\n")

    return xb
end;
