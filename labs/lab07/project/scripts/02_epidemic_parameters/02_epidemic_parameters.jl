using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

N = 17_854.0

I0 = 199.0
R0 = 35.0
S0 = N - I0 - R0

u0 = [S0, I0, R0]

beta = 0.02

alpha_values = [0.5e-5, 1.0e-5, 1.5e-5]

tspan = (0.0, 120.0)
save_step = 0.1

mkpath("data")
mkpath("plots")

results = DataFrame(
    alpha = Float64[],
    peak_infected = Float64[],
    peak_time = Float64[],
    final_susceptible = Float64[],
    final_infected = Float64[],
    final_recovered = Float64[],
)

plot_result = plot(
    xlabel="Time",
    ylabel="Infected",
    title="Influence of infection rate",
)

for alpha in alpha_values

    function epidemic!(du, u, p, t)
        S, I, R = u

        infection = alpha * S * I
        recovery = beta * I

        du[1] = -infection
        du[2] = infection - recovery
        du[3] = recovery
    end

    problem = ODEProblem(
        epidemic!,
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
    susceptible = getindex.(solution.u, 1)
    infected = getindex.(solution.u, 2)
    recovered = getindex.(solution.u, 3)

    peak_index = argmax(infected)

    push!(
        results,
        (
            alpha,
            infected[peak_index],
            time[peak_index],
            last(susceptible),
            last(infected),
            last(recovered),
        ),
    )

    plot!(
        plot_result,
        time,
        infected;
        label="alpha = $(alpha)",
        linewidth=2,
    )
end

CSV.write(
    "data/epidemic_parameter_scan.csv",
    results,
)

savefig(
    plot_result,
    "plots/epidemic_parameter_scan.png",
)

println(results)
println()
println("LAB05 PARAMETERS OK")
