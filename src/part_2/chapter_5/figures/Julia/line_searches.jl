
using ForwardDiff


# Line searches
function uniform(f,a,b,n,k; verbose=false)
    # divides interval [a,b] in n clusters
    # repeats k times around best grid point
    # returns best λ and function evaluation θ
    a_aux = zeros(k+1);
    b_aux = zeros(k+1);
    a_aux[1] = a;
    b_aux[1] = b;

    λ, θ = 0.0, 0.0

    for i = 1:k
        #interval size δ
        δ = (b_aux[i]-a_aux[i])/n
        grid_points = [a_aux[i] + k*δ for k =0:n]
        θ, λ_pos = findmin(f.(grid_points))
        λ = grid_points[λ_pos]
        #new interval
        a_aux[i+1] = λ - δ
        b_aux[i+1] = λ + δ

        if verbose
            println("interval=[", a_aux[i], " ",b_aux[i], "]")
            println("k=",i, ": θ=", θ, "/ λ=", λ)
        end    
    end
    return λ
end;

function dichotomous(f, a, b, l=1e-8, ϵ=1e-10; verbose=false)
    #uses reference points λ < μ spaced by 2e-13.
    #Stops when b - a < l
    #returns minimum point found λ

    bᵢ = b
    aᵢ = a
    while bᵢ - aᵢ > l
        mid = (aᵢ + bᵢ)/2
        λᵢ = mid - ϵ
        μᵢ = mid + ϵ
        θλ = f(λᵢ)
        θμ = f(μᵢ)

        if θλ < θμ
            aᵢ₊₁ = aᵢ
            bᵢ₊₁ = μᵢ
        else
            aᵢ₊₁ = λᵢ
            bᵢ₊₁ = bᵢ
        end
        aᵢ = aᵢ₊₁
        bᵢ = bᵢ₊₁
        
        if verbose
            println("interval=[", aᵢ, " ",bᵢ, "]")
        end    
    end
    return  (aᵢ + bᵢ)/2
end;

function golden_ratio(f, a, b, l=1e-9; verbose = false)
    α = 1/((1+sqrt(5))/2) # φ = golden ratio
    bᵢ, aᵢ = b, a
    aᵢ₊₁, bᵢ₊₁, λᵢ₊₁, μᵢ₊₁, θμᵢ₊₁, θλᵢ₊₁ = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    θλ = f(aᵢ + (1-α)*(bᵢ - aᵢ))
    θμ = f(aᵢ + α*(bᵢ - aᵢ))

    while bᵢ - aᵢ > l
        λᵢ = aᵢ + (1-α)*(bᵢ - aᵢ)
        μᵢ = aᵢ + α*(bᵢ - aᵢ)

        if θλ > θμ # move right
            aᵢ₊₁ = λᵢ
            bᵢ₊₁ = bᵢ
            λᵢ₊₁ = μᵢ
            μᵢ₊₁ = aᵢ₊₁ + α*(bᵢ₊₁ - aᵢ₊₁)
            θμᵢ₊₁ = f(μᵢ₊₁)
            θλᵢ₊₁ = θμ
        else # move left
            aᵢ₊₁ = aᵢ
            bᵢ₊₁ = μᵢ
            μᵢ₊₁ = λᵢ
            λᵢ₊₁ = aᵢ₊₁ + (1-α)*(bᵢ₊₁ - aᵢ₊₁)
            θλᵢ₊₁ = f(λᵢ₊₁)
            θμᵢ₊₁ = θλ
        end
        aᵢ = aᵢ₊₁
        bᵢ = bᵢ₊₁
        θλ = θλᵢ₊₁
        θμ = θμᵢ₊₁
        
        if verbose
            println("interval=[", aᵢ, " ",bᵢ, "]")
        end    
    end
    return  (aᵢ + bᵢ)/2
end;

# Helper function for fibonacci
function find_n(a, b, ϵ=1e-8)
    #returns n such that (b-a)/F_n < 1e-6
    i=1
    previous_F = 0
    current_F = 1
    numbers =[1]
    while (b-a)/(previous_F + current_F) > ϵ
        #calculate Fibonacci numbers without storing
        #or using recursion.
        new_current_F = previous_F + current_F
        previous_F = current_F
        current_F = new_current_F
        push!(numbers, current_F)
        i = i+1
    end
    return i, numbers
end;

function fibonacci(f, a, b, ϵ=1e-8; verbose=false)
    n, seq = find_n(a,b) #quite expensive!
    bᵢ = b
    aᵢ = a
    λᵢ = aᵢ + (seq[n-2]/seq[n])*(bᵢ - aᵢ)
    μᵢ = aᵢ + (seq[n-1]/seq[n])*(bᵢ - aᵢ)
    θλ = f(λᵢ)
    θμ = f(μᵢ)
    i=1
    while i < n-2
        λᵢ = aᵢ + (seq[n-i-2]/seq[n-i])*(bᵢ - aᵢ)
        μᵢ = aᵢ + (seq[n-i-1]/seq[n-i])*(bᵢ - aᵢ)
        if θλ > θμ #move right
            aᵢ₊₁ = λᵢ
            bᵢ₊₁ = bᵢ
            λᵢ₊₁ = μᵢ
            μᵢ₊₁ = aᵢ₊₁ + (seq[n-i-1]/seq[n-i])*(bᵢ₊₁ - aᵢ₊₁)
            θμᵢ₊₁ = f(μᵢ₊₁)
            θλᵢ₊₁ = θμ
        else #move left
            aᵢ₊₁ = aᵢ
            bᵢ₊₁ = μᵢ
            μᵢ₊₁ = λᵢ
            λᵢ₊₁ = aᵢ₊₁ + (seq[n-i-2]/seq[n-i])*(bᵢ₊₁ - aᵢ₊₁)
            θλᵢ₊₁ = f(λᵢ₊₁)
            θμᵢ₊₁ = θλ
        end
        aᵢ = aᵢ₊₁
        bᵢ = bᵢ₊₁
        θλ = θλᵢ₊₁
        θμ = θμᵢ₊₁
        i=i+1
        
        if verbose
            println("interval=[", aᵢ, " ",bᵢ, "]")
        end    
    end
    #last artificial iteration, since ratios become the same (1/2)
    λn = λᵢ
    μn = μᵢ + ϵ
    θλn = f(λn)
    θμn = f(μn)

    if θλn > θμn
        an = λn
        bn = bᵢ
    else
        an = aᵢ
        bn = μn
    end
    return (an + bn)/2
end;

function newton(f, xᵢ, ϵ=1e-8; verbose = false)
    # Precalculates first and second derivative.
    D(f, x)  = ForwardDiff.derivative(f, x)
    # Hack to deal with unidimensional second-order derivative
    D²(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)
    # Newtown step
    while abs(D(f, xᵢ)) > ϵ
        xᵢ₊₁ = xᵢ - D(f, xᵢ)/D²(f, xᵢ)
        xᵢ = xᵢ₊₁
        
        if verbose
            println("Current point: $xᵢ")
        end    
    end
    return xᵢ
end;

function bisection(f, a, b, l=1e-9; verbose = false)
    aᵢ, bᵢ = a, b
    λ, aᵢ₊₁, bᵢ₊₁, ∇f = 0.0, 0.0, 0.0, 0.0
   
    while bᵢ - aᵢ > l
        λ = (aᵢ + bᵢ)/2
        ∇f= ForwardDiff.derivative(f,λ)
        
        if abs(∇f) < l # equivalent to ∇f = 0
            return λ
        elseif ∇f > 0 # left move
            aᵢ₊₁ = aᵢ
            bᵢ₊₁ = λ
        else # move right
            aᵢ₊₁ = λ
            bᵢ₊₁ = bᵢ
        end
        
        aᵢ = aᵢ₊₁
        bᵢ = bᵢ₊₁
        
        if verbose
            println("interval=[", aᵢ, " ",bᵢ, "]")
        end    
    end
        return λ
end;

function armijo(f, α=1e-1, β=0.7, λ=5)

    ∇f0 = ForwardDiff.derivative(f, 0)
    f0 = f(0)
    
    # no function evaluations other than those above!
    while f(λ) > f0 + α*λ*∇f0
        λ = β*λ
    end
    return λ
end;
