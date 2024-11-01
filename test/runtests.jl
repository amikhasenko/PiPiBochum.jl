using PiPiBochum
import PiPiBochum: mpi, c, ρ, amplitude
using Test

@testset "Test c and amplitude functions" begin
    @test imag(c(mpi, mpi, 0.3^2 + 0im)) ≈ -ρ(mpi, mpi, 0.3^2 + 0im)
    @test amplitude(1.1 + 0im) isa Matrix
end
