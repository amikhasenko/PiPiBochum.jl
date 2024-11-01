function Kmatrix(s; pars = default_pars)

    # resonances
    @unpack Gα1, Gα2, Gα3, Gα4, Gα5 = pars
    G = [
        Gα1 * Gα1',
        Gα2 * Gα2',
        Gα3 * Gα3',
        Gα4 * Gα4',
        Gα5 * Gα5',
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

ρ(m1, m2, s) = sqrt(1 - ((m1 + m2)^2 / s)) * sqrt(1 - (m1 - m2)^2 / s)
qa(m1, m2, s) = ρ(m1, m2, s) * sqrt(s) / 2

function c(m1, m2, s)
    qa_val = qa(m1, m2, s)
    term1 = (2 * qa_val / sqrt(s)) * log((m1^2 + m2^2 - s + 2 * sqrt(s) * qa_val) / (2 * m1 * m2))
    term2 = (m1^2 - m2^2) * ((1 / s) - (1 / (m1 + m2)^2)) * log(m1 / m2)
    Σ = (1 / (π)) * (term1 - term2)
    return -Σ
end

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


amplitude_pipi_hat(s; pars = default_pars) =
    amplitude(s; pars)[1, 1] * ρ(mpi, mpi, s)
