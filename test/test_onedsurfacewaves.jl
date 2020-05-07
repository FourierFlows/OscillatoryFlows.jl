module TestOneDSurfaceWaves

using FourierFlows, Test

using OscillatoryFlows.OneDSurfaceWaves

import OscillatoryFlows.OneDSurfaceWaves: updatevars!

devices = (CPU(),)

function instantiate_with_filterkwargs()
    problem = Problem(nx=4, innerK=0.25, outerK=0.5)
    return typeof(problem) <: FourierFlows.Problem
end

function instantiate_velocity_potential()
    problem = Problem(nx=4, innerK=0.25, outerK=0.5)
    potential = VelocityPotential(problem, nz=4)
    return typeof(potential) <: OneDSurfaceWaves.VelocityPotential
end

function calculate_velocity_potential()
    problem = Problem(nx=4, innerK=0.25, outerK=0.5)
    potential = VelocityPotential(problem, nz=4)
    calculate!(potential)
    return typeof(potential) <: OneDSurfaceWaves.VelocityPotential
end

function updatevars!()
    problem = Problem(nx=4, innerK=0.25, outerK=0.5)
    updatevars!(problem)
    return typeof(problem) <: FourierFlows.Problem
end

for dev in devices  
    println("Testing OneDSurfaceWaves on "*string(typeof(dev)))

    @testset "OneDSurfaceWaves" begin
        @test typeof(Problem(nx=4)) <: FourierFlows.Problem
        @test updatevars!()
        @test instantiate_with_filterkwargs()
        @test instantiate_velocity_potential()
        @test calculate_velocity_potential()
    end
end

end # module
