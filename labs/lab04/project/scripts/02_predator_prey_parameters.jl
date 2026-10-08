# # Параметрическое исследование модели «хищник–жертва»
#
# Исследуется влияние коэффициента естественной смертности
# хищников a.
#
# Базовое значение варианта 58: a = 0.38.
#
# Для стационарного состояния:
#
# $$
# x^* = \frac{c}{d},
# \qquad
# y^* = \frac{a}{b}.
# $$

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

# ## Исходные данные

b = 0.043
c = 0.36
d = 0.052

x0 = 6.0
y0 = 23.0

u0 = [x0, y0]

tspan = (0.0, 100.0)
save_step = 0.05

a_values = [0.30, 0.38, 0.46]

mkpath("data")
mkpath("plots")

results = DataFrame(
    a = Float64[],
    stationary_predators = Float64[],
    stationary_prey = Float64[],
    min_predators = Float64[],
    max_predators = Float64[],
    min_prey = Float64[],
    max_prey = Float64[],
)

plot_result = plot(
    xlabel="Time",
    ylabel="Population хищников",
    title="Influence of predator mortality",
)

# ## Параметрический эксперимент

for mortality in a_values
    function model!(du, u, p, t)
        x, y = u
        du[1] = -mortality * x + b * x * y
        du[2] = c * y - d * x * y
    end

    problem = ODEProblem(
        model!,
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

    x_stationary = c / d
    y_stationary = mortality / b

    push!(
        results,
        (
            mortality,
            x_stationary,
            y_stationary,
            minimum(predators),
            maximum(predators),
            minimum(prey),
            maximum(prey),
        ),
    )

    plot!(
        plot_result,
        time,
        predators;
        label="a = $(mortality)",
        linewidth=2,
    )
end

# ## Сохранение результатов

CSV.write(
    "data/predator_prey_parameter_scan.csv",
    results,
)

savefig(
    plot_result,
    "plots/predator_prey_parameter_scan.png",
)

println(results)
println()
println("LAB04 PARAMETERS OK")
