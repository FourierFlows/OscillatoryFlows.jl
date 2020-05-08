push!(LOAD_PATH,"../src/")

using
    Documenter,
    Literate,
    Plots, # to not capture precompilation output
    OscillatoryFlows,
    OscillatoryFlows.OneDWaveEquation
    
#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "OneDWaveEquation/standing_and_propagating_waves.jl",
    "OneDWaveEquation/two_gaussians.jl",
]


for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
    Literate.notebook(example_filepath, OUTPUT_DIR, documenter=true)
end


#####
##### Build and deploy docs
#####

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid Travis CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/OscillatoryFlows.jl/dev/"
)


makedocs(  modules = [OscillatoryFlows],         
           doctest = false,
             clean = true,
         checkdocs = :all,
            format = format,
           authors = "Navid C. Constantinou and Gregory L. Wagner",
          sitename = "OscillatoryFlows.jl",
             pages = Any[
                         "Home"    => "index.md",
                         "Modules" => Any[
                                          "modules/onedwaveequation.md",
                                          "modules/onedsurfacewaves.md",
                                         ],
                         "Examples" => [
                         "OneDWaveEquation" => Any[
                             "generated/standing_and_propagating_waves.md",
                             "generated/two_gaussians.md"
                             ]
                         ],
                        ]
        )

deploydocs(        repo = "github.com/FourierFlows/OscillatoryFlows.jl.git",
               versions = ["stable" => "v^", "v#.#.#"],
           push_preview = true,
           )
