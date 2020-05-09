module TestOneDWaveEquation

using FourierFlows, Test

using OscillatoryFlows: OneDWaveEquation

devices = (CPU(),)

function test_integral_constraints()
    problem = OneDWaveEquation.Problem(; nx=128, Lx=2π, c=1, β=0, dt=0.01)
    
    σ = 1/10 # variance

    ξ₀(x) = exp( - x^2 / (2σ^2) )
    ∫ξ₀ = sqrt(2π*σ^2)

    u₀(x) = 1/2 * exp( - x^2 / (2σ^2) ) + 1/3 * exp( - x^2 / (2σ^2) )
    ∫u₀ = (1/3 + 1/2)*sqrt(2π*σ^2)

    OneDWaveEquation.set_ξ!(problem, ξ₀)
    OneDWaveEquation.set_u!(problem, u₀)

    return problem.vars.integral_constraints.∫u ≈ ∫u₀ && problem.vars.integral_constraints.∫ξ ≈ ∫ξ₀
end

for dev in devices
    println("Testing OneDWaveEquation on "*string(typeof(dev)))

    @testset "OneDWaveEquation" begin
        @test typeof(OneDWaveEquation.Problem(nx=4)) <: FourierFlows.Problem
        @test test_integral_constraints()
    end
end

end # module
