using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

x0 = 100_000.0
y0 = 10_000.0

u0 = [x0, y0]

tspan = (0.0, 1.0)

condition(u, t, integrator) = min(u[1], u[2])

function stop_battle!(integrator)
    terminate!(integrator)
end

cb = ContinuousCallback(condition, stop_battle!)

function regular_forces!(du, u, p, t)
    x, y = u

    du[1] = -0.12x - 0.9y + sin(t)
    du[2] = -0.3x - 0.1y + cos(t)
end

problem1 = ODEProblem(regular_forces!, u0, tspan)

solution1 = solve(
    problem1,
    Tsit5();
    callback=cb,
    abstol=1e-9,
    reltol=1e-9,
    saveat=0.001,
)

t1 = solution1.t
x1 = max.(getindex.(solution1.u, 1), 0.0)
y1 = max.(getindex.(solution1.u, 2), 0.0)

function regular_partisan!(du, u, p, t)
    x, y = u

    du[1] = -0.25x - 0.96y + sin(2t) + 1
    du[2] = -0.25x * y - 0.3y + cos(20t) + 1
end

problem2 = ODEProblem(regular_partisan!, u0, tspan)

solution2 = solve(
    problem2,
    Tsit5();
    callback=cb,
    abstol=1e-9,
    reltol=1e-9,
    saveat=0.00001,
)

t2 = solution2.t
x2 = max.(getindex.(solution2.u, 1), 0.0)
y2 = max.(getindex.(solution2.u, 2), 0.0)

mkpath("data")
mkpath("plots")

df1 = DataFrame(
    time=t1,
    army_x=x1,
    army_y=y1,
)

df2 = DataFrame(
    time=t2,
    army_x=x2,
    army_y=y2,
)

CSV.write("data/lanchester_regular.csv", df1)
CSV.write("data/lanchester_partisan.csv", df2)

p1 = plot(
    t1,
    x1;
    label="Армия X",
    linewidth=2,
    xlabel="Время",
    ylabel="Численность войск",
    title="Регулярные войска",
)

plot!(
    p1,
    t1,
    y1;
    label="Армия Y",
    linewidth=2,
)

savefig(p1, "plots/lanchester_regular.png")

p2 = plot(
    t2,
    x2;
    label="Армия X",
    linewidth=2,
    xlabel="Время",
    ylabel="Численность войск",
    title="Регулярные войска и партизаны",
)

plot!(
    p2,
    t2,
    y2;
    label="Армия Y",
    linewidth=2,
)

savefig(p2, "plots/lanchester_partisan.png")

p3 = plot(
    t1,
    x1 ./ x0;
    label="X — регулярная модель",
    linewidth=2,
    xlabel="Время",
    ylabel="Доля начальной численности",
    title="Сравнение моделей",
)

plot!(
    p3,
    t1,
    y1 ./ y0;
    label="Y — регулярная модель",
    linewidth=2,
)

plot!(
    p3,
    t2,
    x2 ./ x0;
    label="X — модель с партизанами",
    linewidth=2,
    linestyle=:dash,
)

plot!(
    p3,
    t2,
    y2 ./ y0;
    label="Y — модель с партизанами",
    linewidth=2,
    linestyle=:dash,
)

savefig(p3, "plots/lanchester_comparison.png")

println("Регулярная модель:")
println("  время окончания = ", last(t1))
println("  X = ", last(x1))
println("  Y = ", last(y1))

println()

println("Модель с партизанами:")
println("  время окончания = ", last(t2))
println("  X = ", last(x2))
println("  Y = ", last(y2))

println()
println("LAB02 BASE OK")
