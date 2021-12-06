# Examples presented in class - Lecture 1
# written by Fabricio - 07/09/2018

using JuMP, Cbc
# JuMP is for implementing math. programming models;
# Cbc is for solving them.

# Example 1 - resource allocation

# Problem data
i = 1:2 # i=1: Seattle; i=2: San Diego
j = 1:3 # j=1: New York; j=2: Chicago; j=3: Miami

C = [350 600] # Capacities of the factories
D = [325 300 275] # Demand of clients
T = [2.5 1.7 1.8
     3.5 1.9 1.4] # Transportation costs

# Model implementation
# Creates a model and informs the solver to be used
m = Model(solver = CbcSolver())

# Decision variable for the total transported
@variable(m, x[i,j] >= 0)

# Capacity constraint
@constraint(m, cap[i = 1:2], sum(x[i,j] for j = 1:3) <= C[i])
# Demand constraint
@constraint(m, dem[j = 1:3], sum(x[i,j] for i = 1:2) >= D[j])

# Total distribution cost that we want to minimise
@objective(m, Min, sum(T[i,j]*x[i,j] for i=1:2, j=1:3))

println(m) # Prints the mathematical model for debugging
solve(m) # Solve the model

# Prints the optimal solution
println("\nDistribution plan: \n", getvalue(x))

#################################################################






# Example 2 - portfolio optimisation

using Ipopt, Plots
# Ipopt is for solving nonlinear problems - different technology
# Plots is for plotting (!)

# Read the daily prices data from a .csv file.
data = readcsv("prices.csv")
#First row has the names of the stocks
stock_names = data[1, 1:end-2]
#Last two columns has data and US$ rate
prices_data = data[2:end, 1:end-2]

#Having a peek at the data

cs = plot(prices_data, label = stock_names)
plotlyjs()

gui(cs)
#Returns are calculated as (p(t+1) - p(t))/p(t-1)
returns_data = diff(prices_data) ./ prices_data[1:end-1,:]

#Number of days and stocks in data
T, n =  size(returns_data)
#Calculates expected return and covariance
μ, Σ =  mean(returns_data, 1), cov(returns_data)
# Input form the model: minimum average return required.
r_min = 0.2

port = Model(solver=IpoptSolver())

#allocation variables
@variable(port, 0 <= x[1:n] <= 1);

#notice the division by T to correct the average
@constraint(port, sum(μ[j]*x[j] for j=1:n) >= r_min/T);

@objective(port, Min , sum(x[i]*Σ[i,j]*x[j] for i=1:n,j=1:n));

solve(port);

alloc = getvalue(x)
# Should be the same as r_min
ret = μ*alloc*T
# Looking at risk as the st. deviation of the returns.
risk = sqrt(alloc'*Σ*alloc*T)

# Organising the data for plotting. In an array, they are
# considered different series (with different colours)
stock_names = convert(Array{String}, stock_names)
alloc = convert(Array{Float64}, alloc)
# plots a pie chart
pie(stock_names, alloc, legend = false)

#################################################################






# Example 3 - Robust Knapsack Problem

using Distributions, ECOS, LaTeXStrings
# Distributions includes probability distributions,
# ECOS is a solver that we need because of the stucture of the problem
# LatexStrings is to write LaTeX in our plots

N = 18

#Input data
value = [50, 40, 70, 55, 80, 35, 65, 50, 60, 85, 20, 45, 55, 25, 80, 45, 45, 65]
weight_average = [3.5, 4, 5.5, 5, 6, 4.5, 6, 4, 5.5, 7, 5, 4.5, 7, 3.5, 5.5, 3.5, 4, 6.5]
weight_stdev = weight_average*0.3

#Randomly generating initial data using normal distribution
weight_data = Array{Float64}(100,N)
for j= 1:N
   weight_data[:,j] = rand(Normal(weight_average[j], weight_stdev[j]), 100)
end

#Input data: weight limit (capacity)
capacity = 20

# Protection elipsoid: adds 50% of weight as a protection level
P = 0.5*weight_average

# Declaring the model and solving it as a function.
function solve_robust_model(N, value, weight_average, P, Γ)
    m = Model(solver = ECOSSolver()) #creates the model, select the solver

    @variable(m, 0 <= x[1:N] <= 1 ) # creates the binary variables x, one for each item

    @constraint(m, sum(weight_average[j]*x[j] for j = 1:N) + Γ*norm(P'*x) <= capacity) # declare the knapsack constraint

    @objective(m, Max, sum(value[j]*x[j] for j = 1:N)) # declare the objective function

    solve(m)

    return getvalue(x)
end

# This function simulates the item selection againts feasibility.
function feasibility_estimate(solution, repetitions)
    feasible_count = 0
    actual_weight = zeros(18)
    for n= 1:repetitions
        for j= 1:N
           # generate random weights according to distribution
           actual_weight[j] = rand(Normal(weight_average[j], weight_stdev[j]))
        end
        #if total weight more than capacity => problem infeasible.
        if actual_weight'*solution <= capacity
            feasible_count += 1
        end
    end
    return feasible_count/repetitions
end

# generate a range of Gammas to try
Γ_range = linspace(0,1,10)

# storing the results of each run
feas = []
total_value = []
Γ_used = []

# for each Gamma, solve te model and store the results
for Γ in Γ_range
    println(Γ)
    x = solve_robust_model(N, value, weight_average, P, Γ)
    push!(Γ_used, Γ)
    push!(feas, feasibility_estimate(x,5000))
    push!(total_value, value'x)
end

# potting the results from simulation
p1 = plot(Γ_used, feas, xlabel = L"\Gamma", ylabel = "feas. prob.", legend=false)
p2 = plot(Γ_used, total_value, xlabel = L"\Gamma", ylabel = "total value", legend=false, color=:orange)
plot(p1,p2, layout = (2,1))

#################################################################





# Example 4 - Classification

# Generating random multivariate data.
data11 = rand(Normal(8, 4), 100)
data12 = rand(Normal(12, 2), 100)
data1 = [data11 data12]

data01 = rand(Normal(2, 6), 100)
data02 = rand(Normal(5, 2), 100)
data0 = [data01 data02]

scatter(data1[:,1], data1[:,2], label="Pos. obs.")
scatter!(data0[:,1], data0[:,2], label="Neg. obs.")

function classify(data1, data0, δ)
    m = Model(solver = IpoptSolver())

    @variable(m, u[1:100] >= 0)
    @variable(m, v[1:100] >= 0)
    @variable(m, a[1:2])
    @variable(m, b)

    @constraint(m, error1[i=1:100], sum(a[j]*data1[i,j] for j = 1:2) - b - u[i] <= -1)
    @constraint(m, error0[i=1:100], sum(a[j]*data0[i,j] for j = 1:2) - b + v[i] >= 1)

# This parameter trades off accuracy and robustness for the classifier
    @objective(m, Min, sum(u[i] for i=1:100) + sum(v[i] for i=1:100)
        + δ*(sum(a[j]^2 for j=1:2)))

    solve(m);

    return getvalue(a), getvalue(b), getvalue(u), getvalue(v)
end

a, b, u1, v1 = classify(data1, data0, 0.01);

# Calculating the points on plane defined by a and b
x1 = linspace(-10,20,100);
x2 = (b - a[1]*x1)/ a[2];
# boundaries of the slab
x2u = (b + 1 - a[1]*x1)/ a[2];
x2d = (b - 1 - a[1]*x1)/ a[2];
plot!(x1, [x2u x2 x2d],
    xlabel = "x_1",
    ylabel = "x_2",
    label= [ "", "Classifier", ""] ,
    color=[:green :green :green],
    linestyle=[ :dot :solid :dot],
    )

# another classifier, with δ=50
a, b, u2, v2 = classify(data1, data0, 100);
x2 = (b - a[1]*x1)/ a[2];
x2u = (b + 1 - a[1]*x1)/ a[2];
x2d = (b - 1 - a[1]*x1)/ a[2];
plot!(x1, [x2u x2 x2d],
    xlabel = "x_1",
    ylabel = "x_2",
    legend = false,
    color=[:purple :purple :purple],
    linestyle=[:dot :solid :dot],
    )
