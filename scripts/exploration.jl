using PiPiBochum
using LinearAlgebra
using Random
using Plots
using Parameters

theme(:wong2, frame=:box, grid=false)

default_pars = (
    mf=[0.5146109988244556, 0.9062999999986513, 1.23089000002673, 1.461043944511787, 1.696114327468766],
    #
    Gα1=[0.749866997688989, -0.01257099832673861, 0.2753599978535977, -0.1510199937514032, 0.3610299929020451],
    Gα2=[0.06400735441028882, 0.002039993700021009, 0.7741299935288173, 0.5099954460483236, 0.131119996207024],
    Gα3=[-0.2341669602361275, -0.01031664796738707, 0.7228310629513335, 0.1193373160160431, 0.3679219171982366],
    Gα4=[0.01270001206662291, 0.2670000044701449, 0.09214335545338775, 0.02742288751616556, -0.04024795048926635],
    Gα5=[-0.1424226773316178, 0.2277971435654336, 0.1598113086438209, 0.162720778677211, -0.1739657300479793],
    # 
    c00=0.03728069393605827,
    c01=0.0,
    c02=-0.01398000003371962,
    c03=-0.02202999981025169,
    c04=0.01397000015464572,
    c11=0.0,
    c12=0.0,
    c13=0.0,
    c14=0.0,
    c22=0.02349000177968514,
    c23=0.03100999997418123,
    c24=-0.04002991964937379,
    c33=-0.1376928637961125,
    c34=-0.06721849488474475,
    c44=-0.2840099964663654,
    # 
    s0=0.0091125,
    snormAdler=1.0,
    # s0Adler = 0.1139704455925943
)

function Kmatrix(s; pars=default_pars)

    # resonances
    @unpack Gα1, Gα2, Gα3, Gα4, Gα5 = pars
    G = [
        Gα1 * Gα1',
        Gα2 * Gα2',
        Gα3 * Gα3',
        Gα4 * Gα4',
        Gα5 * Gα5'
    ]

    # bare poles  
    @unpack mf = pars
    M = mf

    # couplings
    @unpack c00, c01, c02, c03, c04 = pars
    @unpack c11, c12, c13, c14 = pars
    @unpack c22, c23, c24 = pars
    @unpack c33, c34, c44 = pars
    C = [
        c00 c01 c02 c03 c04;
        c01 c11 c12 c13 c14;
        c02 c12 c22 c23 c24;
        c03 c13 c23 c33 c34;
        c04 c14 c24 c34 c44
    ]

    # Adler term 
    @unpack s0, snormAdler = pars
    adlerTerm = (s - s0) / snormAdler

    # Result
    K = adlerTerm * (G' * (M .^ 2 .- s) .^ (-1) + 5 * C)
    return K
end

#  masses of the decay products 
const mpi = 0.13957
const meta = 0.547862
const m2pi = 0.26996  # m2pi = 2 mpi
const mKp = 0.49367
const mK0 = 0.497614
const metaprime = 0.95778

ρ(m1, m2, s) = sqrt(1 - ((m1 + m2)^2 / s)) * sqrt(1 - (m1 - m2)^2 / s)
qa(m1, m2, s) = ρ(m1, m2, s) * sqrt(s) / 2

function c(m1, m2, s)
    qa_val = qa(m1, m2, s)
    term1 = (2 * qa_val / sqrt(s)) * log((m1^2 + m2^2 - s + 2 * sqrt(s) * qa_val) / (2 * m1 * m2))
    term2 = (m1^2 - m2^2) * ((1 / s) - (1 / (m1 + m2)^2)) * log(m1 / m2)
    Σ = (1 / (π)) * (term1 - term2)
    return -Σ
end

@assert imag(c(mpi, mpi, 0.3^2 + 0im)) ≈ -ρ(mpi, mpi, 0.3^2 + 0im)



ChewMmat(s) = diagm(
    [
    c(mpi, mpi, s),
    c(m2pi, m2pi, s),
    c(mKp, mK0, s),
    c(meta, meta, s),
    c(meta, metaprime, s),
]
)


function amplitude(s; pars=default_pars)
    _K = Kmatrix(s; pars)
    _C = ChewMmat(s)
    𝕀 = Matrix(I, (5, 5))
    inv(𝕀 + _K * _C) * _K
end

@assert amplitude(1.1 + 0im) isa Matrix


let
    ev = range(0.3, 2.0, 1000)
    yv = map(ev) do e
        amplitude(e^2 + 0im)[1, 1]
    end
    plot()
    plot!(ev, yv |> real, lab="real")
    plot!(ev, yv |> imag, lab="imag")
    plot!()
end

amplitude_pipi_hat(s; pars=default_pars) =
    amplitude(s; pars)[1, 1] * ρ(mpi, mpi, s)

let
    ev = range(0.3, 2.0, 5000)
    yv = map(ev) do e
        amplitude_pipi_hat(e^2 + 0im)
    end
    plot()
    plot!(yv, lab="")
    plot!()
end


let
    ev = range(0.3, 2.0, 1000)
    yv = map(ev) do e
        Kmatrix(e^2)[1, 1]
    end
    plot(ev, yv, ylim=(-1, 1))
end
