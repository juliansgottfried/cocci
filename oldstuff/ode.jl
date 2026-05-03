import DifferentialEquations as DE
import Plots

function equations(du, u, p, t)
    S, I1, I2, I12, C1, C2 = u
    r, Kr, β, μ, ν, α, Kc, γ = p

    Nr = S + I1 + I2 + I12
    Nc = C1 + C2

    du[1] = r * Nr * (1 - Nr / Kr) - S * (β * Nc + μ) # dS
    du[2] = β * S * C1 - I1 * (β * C2 + μ) # dI1
    du[3] = β * S * C2 - I2 * (β * C1 + ν) # dI2
    du[4] = β * (I1 * C2 + I2 * C1) - ν * I12 # dI12
    
    du[5] = α * (μ * I1 + ν * I12) * (1 - Nc / Kc) - γ * C1 # dC1
    du[6] = α * ν * (I2 + I12) * (1 - Nc / Kc) - γ * C2 # dC2
end

# S, I1, I2, I12, C1, C2
u0 = [200, 0, 0, 0, 100, 1]

# r, Kr, β, μ, ν, α, Kc, γ
p = [4, 500.0, 0.01, 0.25, 1, 10.0, 1000.0, 25.0]

tspan = (0.0, 50.0)

prob = DE.ODEProblem(equations, u0, tspan, p)
sol = DE.solve(prob, DE.Tsit5(), reltol = 1e-8, abstol = 1e-8)

Plots.plot(sol, linewidth = 1, 
    xaxis = "time",
    yaxis = "size", label = ["S" "I1" "I2" "I12" "C1" "C2"])