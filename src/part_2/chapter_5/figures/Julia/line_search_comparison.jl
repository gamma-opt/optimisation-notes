using Plots, BenchmarkTools

include("line_searches.jl")

# Testing
a,b = -10,10;

# easy function
# f(x) = (x-2)^2 +2;
# hard function
f(x) = log(1+abs(x-4)^(2+sin(x)))

x = range(-10,10,length=1000);
plot(x, f.(x), xaxis=([-10,10], "λ"))

# Comparing the time and solution of each method
print("Uniform: ")
λ1 = uniform(f,a,b,5,5)
println("point found : $(round(λ1, digits=5))")
@btime uniform(f,a,b,n,5)

print("Dichotomous: ")
λ2 = dichotomous(f,a,b)
println("point found : $(round(λ2, digits=5))")
@btime dichotomous(f,a,b)

print("Golden ratio: ")
λ3 = golden_ratio(f,a,b)
println("point found : $(round(λ3, digits=5))")
@btime golden_ratio(f,a,b)

print("Fibonacci: ")
λ4 = fibonacci(f,a,b)
println("point found : $(round(λ4, digits=5))")
@btime fibonacci(f,a,b)

print("Newton: ")
λ5 = newton(f, 1)
println("point found : $(round(λ5, digits=5))")
@btime newton(f, 1)

print("Bisection: ")
λ6 = bisection(f,a,b)
println("point found : $(round(λ6, digits=5))")
@btime bisection(f,a,b)

print("Armijo: ")
λ7 = armijo(f)
println("point found : $(round(λ7, digits=5))")
@btime armijo(f);

