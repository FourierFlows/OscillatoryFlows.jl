using
  FourierFlows,
  Test

# the devices on which tests will run
devices = (CPU(),)
@has_cuda devices = (CPU(), GPU())
@has_cuda using CuArrays

# Run tests
testtime = @elapsed begin

for dev in devices
  
  println("testing on "*string(typeof(dev)))
  
  @testset "empty testset" begin
    @test 1==1
  end
end

end # time

println("Total test time: ", testtime)
