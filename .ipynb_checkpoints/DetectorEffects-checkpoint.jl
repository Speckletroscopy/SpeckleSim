module DetectorEffects

using Random, Distributions, Plots

"""
    adjusted_mean_count_rate(ph::Number, dk::Number, dead::Number)

Returns the adjusted mean count rate.

INPUTS:
    ph: photon count rate
    dk: dark count rate
    dead: detector dead time
"""
function adjusted_mean_count_rate(ph::Number, dc::Number, dead::Number)
    numerator = ph+dc
    denominator = 1+dead*(numerator)
    return numerator/denominator
end

function run()
    dead = collect(1:0.1:10) # 1-10 ns deadtime
    ph = 0.05 # 50 MHz count rate
    dc = 1.0e-6 # 1 KHz dark count rate
    plot(dead,adjusted_mean_count_rate.(ph,dc,dead)/ph)
    
end
export run

end

using .DetectorEffects

DetectorEffects.run()
