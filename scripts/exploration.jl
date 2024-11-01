using SpiralGalaxy
using LinearAlgebra
using Random
using Plots
using Parameters

theme(:wong2, frame = :box, grid = false)

default_pars = (mf = [0.5146109988244556, 0.9062999999986513, 1.23089000002673, 1.461043944511787, 1.696114327468766],
    # 
    gpipif0500 = 0.749866997688989,
    g4pif0500 = -0.01257099832673861,
    gKKf0500 = 0.2753599978535977,
    getaetaf0500 = -0.1510199937514032,
    getaetaprimef0500 = 0.3610299929020451,
    # 
    gpipif0980 = 0.06400735441028882,
    g4pif0980 = 0.002039993700021009,
    gKKf0980 = 0.7741299935288173,
    getaetaf0980 = 0.5099954460483236,
    getaetaprimef0980 = 0.131119996207024,
    # 
    gpipif01370 = -0.2341669602361275,
    g4pif01370 = -0.01031664796738707,
    gKKf01370 = 0.7228310629513335,
    getaetaf01370 = 0.1193373160160431,
    getaetaprimef01370 = 0.3679219171982366,
    # 
    gpipif01500 = 0.01270001206662291,
    g4pif01500 = 0.2670000044701449,
    gKKf01500 = 0.09214335545338775,
    getaetaf01500 = 0.02742288751616556,
    getaetaprimef01500 = -0.04024795048926635,
    # 
    gpipif01710 = -0.1424226773316178,
    g4pif01710 = 0.2277971435654336,
    gKKf01710 = 0.1598113086438209,
    getaetaf01710 = 0.162720778677211,
    getaetaprimef01710 = -0.1739657300479793,
    # 
    c00 = 0.03728069393605827,
    c01 = 0.0,
    c02 = -0.01398000003371962,
    c03 = -0.02202999981025169,
    c04 = 0.01397000015464572,
    c11 = 0.0,
    c12 = 0.0,
    c13 = 0.0,
    c14 = 0.0,
    c22 = 0.02349000177968514,
    c23 = 0.03100999997418123,
    c24 = -0.04002991964937379,
    c33 = -0.1376928637961125,
    c34 = -0.06721849488474475,
    c44 = -0.2840099964663654,
    # 
    s0 = 0.0091125,
    snormAdler = 1.0,
)

function Kmatrix(s; pars = default_pars)

    @unpack gpipif0500, g4pif0500, gKKf0500, getaetaf0500, getaetaprimef0500 = pars
    @unpack gpipif0980, g4pif0980, gKKf0980, getaetaf0980, getaetaprimef0980 = pars
    @unpack gpipif01370, g4pif01370, gKKf01370, getaetaf01370, getaetaprimef01370 = pars
    @unpack gpipif01500, g4pif01500, gKKf01500, getaetaf01500, getaetaprimef01500 = pars
    @unpack gpipif01710, g4pif01710, gKKf01710, getaetaf01710, getaetaprimef01710 = pars

    @unpack mf = pars

    # bare poles  
    M = mf

    # couplings
    G0500 = [gpipif0500, g4pif0500, gKKf0500, getaetaf0500, getaetaprimef0500]
    G0980 = [gpipif0980, g4pif0980, gKKf0980, getaetaf0980, getaetaprimef0980]
    G1370 = [gpipif01370, g4pif01370, gKKf01370, getaetaf01370, getaetaprimef01370]
    G1500 = [gpipif01500, g4pif01500, gKKf01500, getaetaf01500, getaetaprimef01500]
    G1710 = [gpipif01710, g4pif01710, gKKf01710, getaetaf01710, getaetaprimef01710]

    G = [
        G0500 * G0500',
        G0980 * G0980',
        G1370 * G1370',
        G1500 * G1500',
        G1710 * G1710',
    ]

    @unpack c00, c01, c02, c03, c04 = pars
    @unpack c11, c12, c13, c14 = pars
    @unpack c22, c23, c24 = pars
    @unpack c33, c34, c44 = pars
    # 
    C = [
        c00 c01 c02 c03 c04;
        c01 c11 c12 c13 c14;
        c02 c12 c22 c23 c24;
        c03 c13 c23 c33 c34;
        c04 c14 c24 c34 c44
    ]

    # s0Adler = 0.1139704455925943
    @unpack s0, snormAdler = pars
    adlerTerm = (s - s0) / snormAdler

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


function amplitude(s; pars = default_pars)
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
    plot!(ev, yv |> real, lab = "real")
    plot!(ev, yv |> imag, lab = "imag")
    plot!()
end

amplitude_pipi_hat(s; pars = default_pars) =
    amplitude(s; pars)[1, 1] * ρ(mpi, mpi, s)

let
    ev = range(0.3, 2.0, 5000)
    yv = map(ev) do e
        amplitude_pipi_hat(e^2 + 0im)
    end
    plot()
    plot!(yv, lab = "")
    plot!()
end


let
    ev = range(0.3, 2.0, 1000)
    yv = map(ev) do e
        Kmatrix(e^2)[1, 1]
    end
    plot(ev, yv, ylim = (-1, 1))
end
