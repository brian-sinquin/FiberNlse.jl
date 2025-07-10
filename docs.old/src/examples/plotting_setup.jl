using Plots
using FFTW
# Set default plot style
default(
    titlefontcolor=:blue,
    titlefontsize=12,
    titlelocation=:left,
    legendfontsize=10,
    framestyle=:box,
    minorgrid=true
)

# Color schemes
temporal_colors = :viridis
spectral_colors = :magma

"""
    plot_temporal_evolution(sol, tlims=nothing)

Plot temporal evolution of the field amplitude.
"""
function plot_temporal_evolution(sol; title="Temporal Evolution", tlims=nothing)
    p = heatmap(sol.t * 1e12, sol.z * 100,
        abs2.(sol.At) ./ maximum(abs2.(sol.At)),
        xlabel="Time (ps)",
        ylabel="Distance (cm)",
        title=title,
        color=temporal_colors)

    if !isnothing(tlims)
        plot!(p, xlims=tlims)
    end
    return p
end

"""
    plot_spectral_evolution(sol, flims=nothing)

Plot spectral evolution of the field.
"""
function plot_spectral_evolution(sol; title="Spectral Evolution", flims=nothing)
    # Calculate frequency grid
    df = 1 / (sol.t[end] - sol.t[1])
    f = fftshift(fftfreq(length(sol.t), df))

    # Calculate spectral intensity
    Af = abs2.(fftshift(fft(sol.At, 2), 2))
    Af = Af ./ maximum(Af)

    p = heatmap(f * 1e-12, sol.z * 100,
        Af,
        xlabel="Frequency (THz)",
        ylabel="Distance (cm)",
        title=title,
        color=spectral_colors)

    if !isnothing(flims)
        plot!(p, xlims=flims)
    end
    return p
end

"""
    plot_temporal_slices(sol, z_indices, labels)

Plot temporal profiles at specified distances.
"""
function plot_temporal_slices(sol, z_indices, labels; title="Temporal Profiles")
    p = plot(xlabel="Time (ps)",
        ylabel="Normalized Power",
        title=title)

    for (i, idx) in enumerate(z_indices)
        plot!(p, sol.t * 1e12,
            abs2.(sol.At[idx, :]) ./ maximum(abs2.(sol.At[idx, :])),
            label=labels[i])
    end
    return p
end

"""
    plot_spectral_slices(sol, z_indices, labels)

Plot spectral profiles at specified distances.
"""
function plot_spectral_slices(sol, z_indices, labels; title="Spectral Profiles")
    df = 1 / (sol.t[end] - sol.t[1])
    f = fftshift(fftfreq(length(sol.t), df))

    p = plot(xlabel="Frequency (THz)",
        ylabel="Normalized Power",
        title=title)

    for (i, idx) in enumerate(z_indices)
        Af = fftshift(fft(sol.At[idx, :]))
        plot!(p, f * 1e-12,
            abs2.(Af) ./ maximum(abs2.(Af)),
            label=labels[i])
    end
    return p
end

"""
    calculate_fwhm(t, power)

Calculate the Full Width at Half Maximum of a pulse.
"""
function calculate_fwhm(t, power)
    max_power = maximum(power)
    above_half = power .>= max_power / 2
    t_above = t[above_half]
    return (maximum(t_above) - minimum(t_above))
end
