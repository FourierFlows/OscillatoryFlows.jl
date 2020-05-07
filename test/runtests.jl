using FourierFlows, Test

# The devices on which tests will run
devices = (CPU(),)

#@has_cuda devices = (CPU(), GPU())
#@has_cuda using CuArrays

# Run tests
testtime = @elapsed begin
    include("test_onedwaveequation.jl")
    include("test_onedsurfacewaves.jl")
end # time

println("Total test time: ", testtime)
