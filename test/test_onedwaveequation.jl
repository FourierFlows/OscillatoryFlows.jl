module TestOneDWaveEquation

using FourierFlows, Test

using OscillatoryFlows: OneDWaveEquation

devices = (CPU(),)

for dev in devices
    println("Testing OneDWaveEquation on "*string(typeof(dev)))

    @testset "OneDWaveEquation" begin
        @test typeof(OneDWaveEquation.Problem(nx=4)) <: FourierFlows.Problem
    end
end

end # module
