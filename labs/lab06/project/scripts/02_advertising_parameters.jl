# # Параметрическое исследование модели рекламы
#
# Исследуется влияние коэффициента сарафанного радио b
# в модели
#
# dn/dt = (a + b*n)(N - n)

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

N = 761.0
n0 = 6.0

a = 0.00003
b_values = [0.4, 0.82, 1.2]

tspan = (0.0, 0.08)

mkpath("data")
mkpath("plots")

results = DataFrame(
    b = Float64[],
    max_speed = Float64[],
    max_speed_time = Float64[],
    informed_at_max = Float64[],
)

plot_result = plot(
    xlabel="Time",
    ylabel="Informed audience",
    title="Influence of word-of-mouth coefficient",
)

for b in b_values

    function model!(du, u, p, t)
        n = u[1]
        du[1] = (a + b * n) * (N - n)
    end

    problem = ODEProblem(
        model!,
        [n0],
        tspan,
    )

    solution = solve(
        problem,
        Tsit5();
        saveat=0.00001,
        abstol=1e-10,
        reltol=1e-10,
    )

    time = solution.t
    informed = getindex.(solution.u, 1)

    speeds = [
        (a + b * n) * (N - n)
        for n in informed
    ]

    index_max = argmax(speeds)

    push!(
        results,
        (
            b,
            speeds[index_max],
            time[index_max],
            informed[index_max],
        ),
    )

    plot!(
        plot_result,
        time,
        informed;
        label="b = $(b)",
        linewidth=2,
    )
end

CSV.write(
    "data/advertising_parameter_scan.csv",
    results,
)

savefig(
    plot_result,
    "plots/advertising_parameter_scan.png",
)

println(results)
println()
println("LAB06 PARAMETERS OK")
