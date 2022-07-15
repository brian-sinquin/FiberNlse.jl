using Test

tests = [
    "soliton"
]

const testdir = dirname(@__FILE__)


for t in tests
    tp = joinpath(testdir, "$(t).jl")
    @testset "$(t)" begin
      include(tp)
    end
end
