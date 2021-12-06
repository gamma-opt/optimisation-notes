using  Plots, LaTeXStrings
pyplot();

# Include our Newton's method implementation with bisection as line search
include("optimisation_methods.jl")

# We will consider this problem
f(x) = (x[1]-2)^4 + (x[1]-2x[2])^2   # min f(x)
g(x) = x[1]^2 - x[2]                  # g(x) ≦ 0

# Penalty function
B(x) = -1/g(x)
#B(x) = -log(-g(x))


# Initial parameters
μ = 10          # penalty term. Try different values
β = 0.2         # amount of reduction in penalty
ϵ = 1e-7        # tolerance
x₀ = [0.5,1]  # initial feasible point

# Method of multipliers implementation
function barrier_method(μ, f, B, x₀, β, ϵ)
        xᵏ =  x₀
        traj = zeros(2,50)
        traj[:,1] = x₀
        k = 2
        println("\nStarting barrier method...")
        while (abs(μ*B(xᵏ)) > ϵ)
            println(B(xᵏ))
            F(x) = f(x) + μ*B(x)    # barrier problem
            xᵏ,⋅,⋅,⋅ = newton(F,xᵏ) # solver barrier prob. with Newton's
            traj[:,k] = xᵏ          # record trajectory for plotting
            println("Current point: ", xᵏ, "/ Objective : ", f(xᵏ))
            μ = β*μ                # update μ
            k = k+1
        end
        traj = traj[:,1:k-1]
        return traj
end

traj = barrier_method(μ, f, B, x₀, β, ϵ)

# Plotting the contours of the function optimised and the trajectory
n = 5000
x = range(-200,stop=200,length=n);
y = range(-200,stop=200,length=n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

# Plot function contours
contour(x,y,z, # Change title accordingly
        levels = [0.25, 1, 2.5 , 5, 9],
        xaxis = (L"$x_1$", (0,2)),
        yaxis = (L"$x_2$", (0,2)),
        clims = (0,12),
        #clabels = true,
        legend = false,
        aspect_ratio = :equal)

# Plot feasible region
plot!(x, x.^2,
        #fill= (10,0.1),
        color = 1)
annotate!([(0.25, 1.75, text(L"$x_1^2 - x_2 \leq 0$",9,:left,
                RGBA(0.0,0.605603,0.97868,1.0)))])

# Plot trajectory

plot!(traj[1,:], traj[2,:], marker=:o)

annotate!([(0.5, 0.9, text(L"$x_0$",9,:left,
                RGBA(0.2422242978521988,0.6432750931576304,0.30444865153411527,1.0) ))])

# gui()
savefig("ex3-barrier.pdf")
