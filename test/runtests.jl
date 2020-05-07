using
  FourierFlows,
  Test

using OscillatoryFlows: OneDWaveEquation, OneDSurfaceWaves

# the devices on which tests will run
devices = (CPU(),)
#@has_cuda devices = (CPU(), GPU())
#@has_cuda using CuArrays

# Run tests
testtime = @elapsed begin

for dev in devices
  
  println("testing on "*string(typeof(dev)))
 
  @testset "OneDWaveEquation" begin
      @test typeof(OneDWaveEquation.Problem(nx=4)) <: FourierFlows.Problem
  end

  @testset "OneDSurfaceWaves" begin
      @test typeof(OneDSurfaceWaves.Problem(nx=4)) <: FourierFlows.Problem
  end
end

end # time

println("Total test time: ", testtime)
