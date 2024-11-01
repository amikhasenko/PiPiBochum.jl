using PiPiBochum
using Plots
using Parameters

import PiPiBochum: amplitude, amplitude_pipi_hat
import PiPiBochum: default_pars_paper


let
    ev = range(0.3, 2.0, 5000)
    yv = map(ev) do e
        amplitude_pipi_hat(e^2 + 0im)
    end
    plot()
    plot!(yv, lab = "")
    yv = map(ev) do e
        amplitude_pipi_hat(e^2 + 0im; pars = default_pars_paper)
    end
    plot!(yv, lab = "")
    plot!()
end

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