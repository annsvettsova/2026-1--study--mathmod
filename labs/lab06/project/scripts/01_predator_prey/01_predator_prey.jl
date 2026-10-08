using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

a = 0.38
b = 0.043
c = 0.36
d = 0.052

x0 = 6.0
y0 = 23.0

u0 = [x0, y0]

tspan = (0.0, 100.0)
save_step = 0.05

mkpath("data")
mkpath("plots")

x_stationary = c / d
y_stationary = a / b

println("Stationary state:")
println("x* = ", x_stationary)
println("y* = ", y_stationary)

function predator_prey!(du, u, p, t)
    x, y = u

    du[1] = -a * x + b * x * y
    du[2] =  c * y - d * x * y
end

problem = ODEProblem(
    predator_prey!,
    u0,
    tspan,
)

solution = solve(
    problem,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

time = solution.t
predators = getindex.(solution.u, 1)
prey = getindex.(solution.u, 2)

results = DataFrame(
    time=time,
    predators=predators,
    prey=prey,
)

CSV.write(
    "data/predator_prey.csv",
    results,
)

p1 = plot(
    time,
    predators;
    label="Predators x(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Population",
    title="Predator-prey model",
)

plot!(
    p1,
    time,
    prey;
    label="Prey y(t)",
    linewidth=2,
)

savefig(
    p1,
    "plots/predator_prey_time.png",
)

p2 = plot(
    prey,
    predators;
    label="Phase trajectory",
    linewidth=2,
    xlabel="Population жертв y",
    ylabel="Population хищников x",
    title="Phase portrait",
)

scatter!(
    p2,
    [y_stationary],
    [x_stationary];
    label="Stationary point",
    markersize=6,
)

savefig(
    p2,
    "plots/predator_prey_phase.png",
)

p3 = plot(
    time,
    predators;
    label="x(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Predators",
    title="Population хищников",
)

savefig(
    p3,
    "plots/predators.png",
)

p4 = plot(
    time,
    prey;
    label="y(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Prey",
    title="Population жертв",
)

savefig(
    p4,
    "plots/prey.png",
)

println()
println("Minimum predators = ", minimum(predators))
println("Maximum predators = ", maximum(predators))
println("Minimum prey = ", minimum(prey))
println("Maximum prey = ", maximum(prey))

println()
println("LAB04 BASE OK")
