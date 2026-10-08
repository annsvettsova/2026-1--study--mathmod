# # Лабораторная работа № 6
#
# ## Модель распространения рекламы
#
# Вариант 58.
#
# Объём потенциальной аудитории:
#
# $$
# N = 761.
# $$
#
# В начальный момент о товаре знают:
#
# $$
# n(0)=6.
# $$

using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

# ## Initial data

N = 761.0
n0 = 6.0

mkpath("data")
mkpath("plots")

# ## Case 1
#
# $$
# \frac{dn}{dt}
# =
# (0.82+0.00003n)(N-n).
# $$

function advertising1!(du, u, p, t)
    n = u[1]

    du[1] = (0.82 + 0.00003 * n) * (N - n)
end

problem1 = ODEProblem(
    advertising1!,
    [n0],
    (0.0, 10.0),
)

solution1 = solve(
    problem1,
    Tsit5();
    saveat=0.01,
    abstol=1e-10,
    reltol=1e-10,
)

t1 = solution1.t
n1 = getindex.(solution1.u, 1)

# ## Case 2
#
# $$
# \frac{dn}{dt}
# =
# (0.00003+0.82n)(N-n).
# $$

function advertising2!(du, u, p, t)
    n = u[1]

    du[1] = (0.00003 + 0.82 * n) * (N - n)
end

problem2 = ODEProblem(
    advertising2!,
    [n0],
    (0.0, 0.05),
)

solution2 = solve(
    problem2,
    Tsit5();
    saveat=0.00001,
    abstol=1e-10,
    reltol=1e-10,
)

t2 = solution2.t
n2 = getindex.(solution2.u, 1)

# ## Maximum advertising speed for case 2
#
# Let
#
# $$
# v(n)=(a+bn)(N-n),
# $$
#
# where
#
# $$
# a=0.00003,\qquad b=0.82.
# $$
#
# Maximum speed occurs at
#
# $$
# n_*=\frac{bN-a}{2b}.
# $$

a2 = 0.00003
b2 = 0.82

n_speed_max = (b2 * N - a2) / (2 * b2)

K = a2 + b2 * N
y0 = a2 + b2 * n0

t_speed_max = log((K - y0) / y0) / K

speed_max =
    (a2 + b2 * n_speed_max) *
    (N - n_speed_max)

println("Case 2 maximum speed:")
println("n* = ", n_speed_max)
println("t* = ", t_speed_max)
println("maximum dn/dt = ", speed_max)

# Dense speed plot around the maximum

speed_time = range(
    0.0,
    0.03;
    length=2000,
)

speed_n = [
    solution2(t)[1]
    for t in speed_time
]

speed_values = [
    (a2 + b2 * n) * (N - n)
    for n in speed_n
]

# ## Case 3
#
# $$
# \frac{dn}{dt}
# =
# (0.2\sin t+0.8\cos(t)n)(N-n).
# $$

function advertising3!(du, u, p, t)
    n = u[1]

    du[1] =
        (0.2 * sin(t) + 0.8 * cos(t) * n) *
        (N - n)
end

problem3 = ODEProblem(
    advertising3!,
    [n0],
    (0.0, 1.5),
)

solution3 = solve(
    problem3,
    Tsit5();
    saveat=0.001,
    abstol=1e-10,
    reltol=1e-10,
)

t3 = solution3.t
n3 = getindex.(solution3.u, 1)

# ## Save data

CSV.write(
    "data/advertising_case1.csv",
    DataFrame(
        time=t1,
        informed=n1,
    ),
)

CSV.write(
    "data/advertising_case2.csv",
    DataFrame(
        time=t2,
        informed=n2,
    ),
)

CSV.write(
    "data/advertising_case3.csv",
    DataFrame(
        time=t3,
        informed=n3,
    ),
)

CSV.write(
    "data/advertising_case2_speed.csv",
    DataFrame(
        time=collect(speed_time),
        informed=speed_n,
        speed=speed_values,
    ),
)

# ## Case 1 plot

p1 = plot(
    t1,
    n1;
    label="n(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Informed audience",
    title="Advertising model: case 1",
)

savefig(
    p1,
    "plots/advertising_case1.png",
)

# ## Case 2 plot

p2 = plot(
    t2,
    n2;
    label="n(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Informed audience",
    title="Advertising model: case 2",
)

savefig(
    p2,
    "plots/advertising_case2.png",
)

# ## Case 2 speed plot

p_speed = plot(
    speed_time,
    speed_values;
    label="dn/dt",
    linewidth=2,
    xlabel="Time",
    ylabel="Advertising speed",
    title="Case 2: advertising speed",
)

scatter!(
    p_speed,
    [t_speed_max],
    [speed_max];
    label="Maximum",
    markersize=6,
)

savefig(
    p_speed,
    "plots/advertising_case2_speed.png",
)

# ## Case 3 plot

p3 = plot(
    t3,
    n3;
    label="n(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Informed audience",
    title="Advertising model: case 3",
)

savefig(
    p3,
    "plots/advertising_case3.png",
)

println()
println("N = ", N)
println("n(0) = ", n0)
println("Case 1 final n = ", last(n1))
println("Case 2 final n = ", last(n2))
println("Case 3 final n = ", last(n3))

println()
println("LAB06 BASE OK")
