using Test
using DSP
using FiberNlse

const TEST_DIR = dirname(@__FILE__)

bundles = ["SSFM"]

for b in bundles 
    @info "Testing $(b)  bundle"
    tests =  readdir(joinpath(dirname(@__FILE__),b), join=true)
    for t in tests
        include(t)
    end
end
