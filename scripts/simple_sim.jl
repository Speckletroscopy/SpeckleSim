### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f1c9d31f-c565-4e96-8801-22352b014397
import Pkg

# ╔═╡ f6daaa32-fe10-4a20-b30e-ba31b63c6b90
Pkg.add("Distributions")

# ╔═╡ 338200b2-c10c-4e7d-8a24-bfd32e6757e5
Pkg.add("DataFrames")

# ╔═╡ 525ba450-c635-4bb9-a2f3-aeb12d5c30d3
Pkg.add("CSV")

# ╔═╡ 64adb0f4-b695-4be1-9b02-c42562951b39
Pkg.add("Plots")

# ╔═╡ 55405eea-1c12-4060-81f5-cbcc082a992c
Pkg.add("IterTools")

# ╔═╡ 7e392b1a-31a8-4364-ac93-742221501f1c
Pkg.add("PlutoUI")

# ╔═╡ b708f1d2-e8dd-40ee-92fd-cb4ae02a6a9a
Pkg.add("StatsBase")

# ╔═╡ 34a840a1-7050-4719-a81e-c4545551ec7d
using Random, Distributions, StatsBase

# ╔═╡ 83920316-034e-4ca1-94fa-5aea4d4d117b
using LinearAlgebra

# ╔═╡ 8287cc85-37d3-4b53-ab65-0aef5d996a72
using Plots

# ╔═╡ 1196f179-17ff-42d9-a582-1b8498b7c372
using IterTools

# ╔═╡ 2d8d5528-d64d-4026-b5a6-93b66e17e609
using PlutoUI

# ╔═╡ c3929637-0224-4f95-a8ce-b1a159084014
using DataFrames

# ╔═╡ 2f28edec-8461-4eda-a95f-05f4db5b3a81
using CSV

# ╔═╡ 0a435f6c-9f0a-11eb-114f-fb8c733dafbb
md"# Brute Force Simulation"

# ╔═╡ 5c181bf0-f92e-40ac-9193-62c48927f60f
md"""
## Simulation Configuration
"""

# ╔═╡ 0fe833e8-a359-41e2-ba4a-b8313a222182
md"""
### Classical field parameters

$(@bind bigN Slider(1:1000, show_value = true))

Number of emitting atoms


$(@bind temp Slider(293:3000, show_value = true))

Temperature (K)

$(@bind tres NumberField(10:500))

Time resolution (ps)

$(@bind tmax NumberField(1.0:0.1:10.0))

Maximum time (ns)

"""

# ╔═╡ c0e353c3-53e7-4419-b8cf-878f6fbc761d
md"""
### Photon count parameters

$(@bind nbar NumberField(1.0:1000.0;default=10))

Expected counts in time tmax

$(@bind ntot NumberField(1:1000000;default=1000))

Total photon counts overall (approximate)
"""

# ╔═╡ 402b4a2a-d22f-4ee7-ae3d-ca5fe7fda110
md"""
Number of classical simulations: $(convert(Integer,ceil(ntot/nbar)))
"""

# ╔═╡ 89b17676-6496-4d37-929d-3f626e329cdd
md"""
### Plotting options

**Available Plots:**	
- Single instance of classical calculations $(@bind classicalPlots CheckBox())	
- Average of all classical calculations $(@bind classicalAvgPlots CheckBox())
- Total counts from all time windows $(@bind countsPlot CheckBox())
- Autocorrelation of photon count time series  $(@bind corrPlot CheckBox())

#### Make Plots [$(@bind makePlots CheckBox())]

"""

# ╔═╡ 30f9f516-568b-4e73-8aa2-4f5bb48f9ab4
md"""
## Plots
"""

# ╔═╡ 57da1ea1-5a5d-49d1-88f2-2624a14f8dff
if makePlots && classicalPlots
	md"""
	### Single instance of classical calculations:
	"""
end

# ╔═╡ b2def348-d4c0-4461-bee6-0373673a8fd8
md"**Download data for classical field and correlation plots:**"

# ╔═╡ 1e48204e-23d8-4581-a160-7ba84b83b8d8
if makePlots && classicalAvgPlots
	md"""
	### Average of all classical calculations:
	"""
end

# ╔═╡ 0f36d9b7-6ba6-402d-9601-53b16e952110
if makePlots && countsPlot
	md"""
	### Total counts:
	"""
end

# ╔═╡ c54073d5-b58a-46ed-b086-d6347c1dbcda
if makePlots && corrPlot
	md"""
	### Autocorrelation of photon counts
	
	Note that the τ = 0 bin matches the number of simulations. This only happens because we are using the autocorrelation for a single beam. This will not be true when I add a beam splitter.
	"""
end

# ╔═╡ 13e70af4-146d-43e7-9fd6-ec8f901185d4
md"""
## Calculations
"""

# ╔═╡ 807e056a-3c89-46f9-8dc2-5dde04023ce2
md"""
### Classical calculations
"""

# ╔═╡ f2e2b9fd-9d21-410c-89d8-2286b163bb2a
# these parameters are needed for all simulations
begin
	# Balmer-α lines specified here
	ωM = 2*π*[456811.0, 456812.0]
	# magnitude of each line
	mag = convert(Vector{ComplexF64},ones(length(ωM)))
	# calculate line differences
	ΔM = ωM .- ωM[1]
	# generate times in tres ps intervals up to 2*tmax	
	times = collect(0:tres*1e-3:2*tmax); 
	# limit the window to tmax to avoid correlation cutoff
	window = convert(Integer,floor(length(times)/2));
	# τ is just the times up to our window
	τ = times[1:window];
end;

# ╔═╡ 7ac53f6b-15ee-440d-acae-005dc3acce51
md"""
f0 (GHz) = $(ωM[1]/(2*π))

Δf (GHz) = $(round((ωM[2]-ωM[1])/(2*π),digits=2))

σ (GHz) = $(round(sqrt(9.178e-14*temp)*ωM[1]/(2*π),digits=2))
"""

# ╔═╡ a3cd4171-5bc5-4988-a5ef-428dd6fa31ab
if makePlots
	md"""
	Size: $(@bind windowSize Slider(1:window,default=window))
	"""
end

# ╔═╡ e2fab25a-ccde-44e8-a3bf-e402fc32d2ba
if makePlots
	md"""
	Position: $(@bind windowPos Slider(0:(window-windowSize)))
	"""
end

# ╔═╡ 8b607526-47b0-4a8a-aa16-2aee351b590f
if makePlots
	md"""
	**Window:**

	Window size (ns) = $(tres*windowSize/1000)

	Window position (ns) = $(windowPos*tres/1000)
	"""
end

# ╔═╡ 21e16aae-1b50-4f31-936a-17525352afbe
if makePlots
	# our calculated g2τ
	g2τCalc = 1
	Emag2 = real.(mag .* conj(mag))
	Emag4 = Emag2 .* Emag2
	sumEmag2 = sum(Emag2)
	sumEmag4 = sum(Emag4)
	term2 = sumEmag4/(bigN*sumEmag2^2)
	g2τCalc -= term2
	
	term3 = sum(Emag2 .* exp.(-im*ΔM .* transpose(τ)),dims=1)/sumEmag2
	term3 = real.(term3 .* conj(term3))
	
	kbOverMhC2 = 9.178e-14;
	σ = sqrt(kbOverMhC2*temp)*ωM[1]
	stauAvg = transpose(bigN .+ bigN*(bigN-1)*exp.(-σ^2*τ .^2))
	term3 = term3 .* stauAvg/bigN^2
	
	g2τCalc =  g2τCalc .+ term3
end;

# ╔═╡ 173ca217-29e6-4c5f-8995-0bd5b62c32dd
md"""
- Our calculated version

$g^{(2)}(\tau) =1-\frac{1}{N}\frac{\sum_{m=1}^M|\mathcal{E}|_m^4}{\left(\sum_{m=1}^M|\mathcal{E}_m|^2\right)^2}+\left|\frac{\sum_{m=1}^M|\mathcal{E}_m|^2e^{-i\Delta_m\tau}}{\sum_{m=1}^M|\mathcal{E}|_m^2}\right|^2\frac{\langle S(\tau)\rangle_\omega}{N^2}$
"""

# ╔═╡ 689bc5a6-2469-41eb-99ce-c300dd73874e
md"""
### Calculate photon counts
"""

# ╔═╡ f2a04b1f-f61c-45d4-a67a-1609f53d452e
md"""
Photon counts are calculated by treating the intensity in each time bin as the average photon count rate, then sampling from a poisson distribution with that average count rate.
"""

# ╔═╡ 730e47f5-1c6a-4fd3-a83f-f877b14cde06
md"""
### Calculate autocorrelation of photon counts

Note that I only look at whether two bins both have counts or not when calculating the autocorrelation. I **do not** look at *how many* counts there are in each bin.
"""

# ╔═╡ 282a7d64-8d67-4a10-9aa9-d64be7c5dd54
md"""
## Functions
"""

# ╔═╡ d2ec427e-50f8-47d2-9a44-cb0d5c0dcb9b
"""
	function countDeltaTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function countTimes(times::Vector,γCounts::Vector)
	out = Vector{Real}(undef,0)
	for (i,counts) in enumerate(γCounts)
		if counts != 0
			countTimes = times[i]*ones(counts)
			out = vcat(out,countTimes)
		end
	end
	return out
end



# ╔═╡ 26e08020-6775-431f-b33a-0efaad94189a
"""
	eFieldParams64(mag::Array{S},ΔM::Array{T},ωN::Array{T},ϕ::Array{S} ) where {T<:Real, S<:Complex}

Static parameters for the electric field
"""
struct eFieldParams
	mag::Vector
	ΔM::Vector
	ωN::Vector
	ϕ::Matrix
	
	function eFieldParams(mag::Vector,ΔM::Vector,ωN::Vector,ϕ::Matrix )
		@assert length(mag) == length(ΔM) "Number of magnitudes must match number of emission lines"
		@assert (size(ϕ)[1] == length(ΔM) && size(ϕ)[2] == length(ωN)) "Must have a unique phase for each n and m"
		new(
			convert(Vector{Complex},mag),
			convert(Vector{Real},ΔM),
			convert(Vector{Real},ωN),
			convert(Matrix{Real},ϕ)
		)
		
	end
end

# ╔═╡ 8175ac47-b261-4942-a70d-f127f22aee74
"""
	function electricField(t::Real,params::eFieldParams)

Returns the electric field value at time t
"""
function electricField(t::Real,params::eFieldParams)
	# generate frequencies
	ωNM = transpose(params.ωN) .+ params.ΔM
	# add the phase
	exponentnm = -im*(t*ωNM+params.ϕ)
	# put them in the exponent
	enm = exp.(exponentnm)
	# multiply by the field magnitude
	fieldnm = params.mag .* enm
	# add it all together
	return sum(ivec(fieldnm))
end

# ╔═╡ d8c63dd0-d6dd-4234-8515-caf5803af683
"""
	function ωnDoppler(ω0::Real,N::Integer,temp::Real,seed::Integer = -1)

Generates N doppler shifted frequencies around frequency ω0 for a source at temperature temp. Seed optional for reproducible results.
"""
function ωnDoppler(ω0::Real,N::Integer,temp::Real,seed::Integer = -1)
	rng = MersenneTwister()
	if seed != -1
		rng = MersenneTwister(seed)
	end
	kbOverMhC2 = 9.178e-14;
	σ = sqrt(kbOverMhC2*temp)*ω0
	d = Normal(ω0,σ)
	return rand(rng,d,N)
end

# ╔═╡ b237ca1e-1719-420c-820d-67b8e18c7a57
"""
	function correlate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}

Calculates correlation between vectors u and v with given offset. Specify averaging window to limit range of correlation. If the window extends beyond the end of one vector, it treats out-of-bounds indices as zero.
"""
function correlate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}
	
	@assert offset <= length(u) "Offset out of bounds"
	@assert window <= length(u) && window <= length(v) "Window must be smaller than input vector lengths"
	
	if window == -1
		window = length(u)
	end
	
	v1 = view(u,1:window)
	v2 = view(v,1+offset:min(window+offset,length(v)))
	if window+offset > length(v)
		v2 = vcat(v2,zeros(window+offset-length(v)))
	end
	
	return dot(v1,v2)/window
end

# ╔═╡ 3620f708-3eac-4cd8-90db-72ac383b5917
"""
	function autocorrelate(u::Vector{T},offset::Integer, window::Integer = -1) where {T<:Number}

Calculates correlation of vector u with itself.
"""
function autocorrelate(u::Vector{T},offset::Integer, window::Integer = -1) where {T<:Number}
	correlate(u,u,offset,window)
end

# ╔═╡ 18f7bc41-751f-4283-8a90-af3ad217b58b
# if ONLY the classical plots are desired, then this calculates a single instance of the classical fields and correlations
if makePlots && (classicalPlots && !(classicalAvgPlots || countsPlot || corrPlot))
	# generate N doppler broadened frequencies
	ω1Doppler = ωnDoppler(ωM[1],bigN,temp)
	# generate N*M random phases
	ϕ1nm = 2*π*rand(Float64,(length(ωM),bigN))
	# construct field parameter object
	testParams1 = eFieldParams(mag,ΔM,ω1Doppler,ϕ1nm)
	
	# calculate electric field vs time
	e1fieldt = map(t->electricField(t,testParams1),times)
	
	# calculate g1τ
	g1τ1Norm = correlate(e1fieldt,conj.(e1fieldt),0,window);
	g1τ1 = abs.(map(i->correlate(e1fieldt,conj.(e1fieldt),i,window),collect(0:window-1))/g1τ1Norm)
	
	# calculate the intensity vs time
	intensity1t = real.( e1fieldt .* conj(e1fieldt))
	
	# calculate g2τ
	g2τ1Norm = mean(intensity1t)^2
	g2τ1 = map(i->autocorrelate(intensity1t,i,window),collect(0:window-1))/g2τ1Norm	
end;

# ╔═╡ 39e8d2c6-c8de-4e0a-b531-60dbaadd6aa1
# the classical average plots and photon-based plots require multiple instances of the classical field calculations 
if makePlots && (classicalAvgPlots || countsPlot || corrPlot)
	# calculate number of trials from the total desired photon counts and the average photon count per trial
	nTrials = convert(Integer,ceil(ntot/nbar))
	# make an array of the average photon counts per trial for array broadcasting
	nPerTrial = bigN*ones(Integer,nTrials)
	# calculate nTrials instances of doppler broadened frequencies
	ωDoppler = ωnDoppler.(ωM[1],nPerTrial,temp);
	# calculate nTrials instances of random phases
	ϕnm = map(n->2*π*rand(Float64,(length(ωM),n)),nPerTrial);
	# generate nTrials instances of field parameters
	testParams = map((x,y)->eFieldParams(mag,ΔM,x,y),ωDoppler,ϕnm);
	
	# calculate nTrials instances of the time dependent electric field
	efieldt = map(x->map(t->electricField(t,x),times),testParams);
	
	# calculate nTrials instances of g1τ
	g1τNorm = map(eft->correlate(eft,conj.(eft),0,window),efieldt);
	g1τ = map((eft,normG1τ)->abs.(map(i->correlate(eft,conj.(eft),i,window),collect(0:window-1))/normG1τ),efieldt,g1τNorm);
	
	# calculate nTrials instances of the time dependent intensity
	intensityt = map(eft->real.( eft .* conj(eft)),efieldt);
	
	# calculate nTrials instances of g2τ
	g2τNorm = map(intens->mean(intens)^2,intensityt)
	g2τ = map((intens,normG2τ)->map(i->autocorrelate(intens,i,window),collect(0:window-1))/normG2τ,intensityt,g2τNorm)
	
	# pick out one instance of everything for plotting
	e1fieldtM = efieldt[1]
	g1τ1M = g1τ[1]
	intensity1tM = intensityt[1]
	g2τ1M = g2τ[1]
end;

# ╔═╡ 9e013576-7682-4d6c-8caa-e3f2c1034f73
if makePlots && classicalPlots
	if !(classicalAvgPlots || countsPlot || corrPlot)
		e1fieldPlot = plot(times[1+windowPos:windowSize+windowPos],real.(e1fieldt)[1+windowPos:windowSize+windowPos],label=false,title="E(t)")
		xlabel!(e1fieldPlot,"t (ns)")

		eCorr1Plot = plot(τ[1+windowPos:windowSize+windowPos],g1τ1[1+windowPos:windowSize+windowPos],label=false,title="g1(τ)")
		xlabel!(eCorr1Plot,"τ (ns)")

		intensity1tPlot = plot(times[1+windowPos:windowSize+windowPos],intensity1t[1+windowPos:windowSize+windowPos],label=false,color=2,title="I(t)")
		xlabel!(intensity1tPlot,"t (ns)")

		iCorr1Plot = plot(τ[1+windowPos:windowSize+windowPos],g2τ1[1+windowPos:windowSize+windowPos],label="data",color=2,title="g2(τ)")
		plot!(iCorr1Plot,τ[1+windowPos:windowSize+windowPos],transpose(g2τCalc)[1+windowPos:windowSize+windowPos],label="calc")
		xlabel!(iCorr1Plot,"τ (ns)")

		plot(e1fieldPlot,intensity1tPlot,eCorr1Plot,iCorr1Plot,layout=(2,2),size=(800,800))
	else
		e1fieldPlot = plot(times[1+windowPos:windowSize+windowPos],real.(e1fieldtM)[1+windowPos:windowSize+windowPos],label=false,title="E(t)")
		xlabel!(e1fieldPlot,"t (ns)")

		eCorr1Plot = plot(τ[1+windowPos:windowSize+windowPos],g1τ1M[1+windowPos:windowSize+windowPos],label=false,title="g1(τ)")
		xlabel!(eCorr1Plot,"τ (ns)")

		intensity1tPlot = plot(times[1+windowPos:windowSize+windowPos],intensity1tM[1+windowPos:windowSize+windowPos],label=false,color=2,title="I(t)")
		xlabel!(intensity1tPlot,"t (ns)")

		iCorr1Plot = plot(τ[1+windowPos:windowSize+windowPos],g2τ1M[1+windowPos:windowSize+windowPos],label="data",color=2,title="g2(τ)")
		plot!(iCorr1Plot,τ[1+windowPos:windowSize+windowPos],transpose(g2τCalc)[1+windowPos:windowSize+windowPos],label="calc",color=3)
		xlabel!(iCorr1Plot,"τ (ns)")

		plot(e1fieldPlot,intensity1tPlot,eCorr1Plot,iCorr1Plot,layout=(2,2),size=(800,800))
	end
end

# ╔═╡ b5342349-5b28-4918-bedb-effa46df10c3
if makePlots && classicalPlots
	if !(classicalAvgPlots || countsPlot || corrPlot)
		classicalFieldDf = DataFrame(
			time              = times, 
			realElectricField = real.(e1fieldt),
			imagElectricField = imag.(e1fieldt),
			intensity         = intensity1t
			)
		classicalCorrDf = DataFrame(
			tau               = τ,
			g1tau             = g1τ1,
			g2tau             = g2τ1			
			)
	else
		classicalFieldDf = DataFrame(
			time              = times, 
			realElectricField = real.(e1fieldtM),
			imagElectricField = imag.(e1fieldtM),
			intensity         = intensity1tM
			)
		classicalCorrDf = DataFrame(
			tau               = τ,
			g1tau             = g1τ1M,
			g2tau             = g2τ1M
			)
	end
	classicalFieldString = string(collect(CSV.RowWriter(classicalFieldDf))...)
	classicalCorrString  = string(collect(CSV.RowWriter(classicalCorrDf))...)
	md"""
	E-field/Intensity:
	$(DownloadButton(classicalFieldString,"classical_field_single.csv"))
	
	Correlations:
	$(DownloadButton(classicalCorrString,"classical_correlations_single.csv"))
	"""
end

# ╔═╡ a3e82734-c1d3-4e9b-8603-ce6bc1cb5018
"""
	function singleDeltaTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function singleDeltaTimes(τ::Vector,γCounts::Vector)
	return τ[map(i->autocorrelate(γCounts,i,length(τ)) > 0 ? true : false,collect(0:length(τ)-1)) ]
end


# ╔═╡ 48fd09eb-f5a0-4bc4-bdbf-3648366978a1
"""
	function γIntensity(intensity::Vector,nbar::Real)

Calculates the photon count rate in each bin of an intensity histogram
"""
function γIntensity(intensity::Vector,nbar::Real)
	nintensity = intensity/sum(intensity)
	return nbar*nintensity
end

# ╔═╡ 31b92f9f-6be0-48cb-9a7e-bfc5edb79bd2
"""
	function poissonCount(nbar::Real)

Returns Poisson distributed counts for average count rate nbar
"""
function poissonCount(nbar::Real)
	p = exp(-nbar)
	s = p
	r = rand()
	count = 0
	while r > s
		count += 1
		p *= nbar/count
		s += p
	end
	return count
end

# ╔═╡ 045fbba2-5ae7-4a41-9e8e-343a5ab41c7b
if makePlots && (countsPlot || corrPlot)
	γCounts = map(intens->poissonCount.(γIntensity(intens,nbar*2)),intensityt)
	γcountTimes = map(γCt->countTimes(times,γCt),γCounts)
	flatγCountTimes = vcat(γcountTimes...)
end;

# ╔═╡ 92d93e89-1d36-4bf0-b3f8-d4370e08857b
if makePlots && countsPlot
	countHist =  fit(Histogram,flatγCountTimes, vcat(τ,tres/10)[1+windowPos:windowSize+windowPos],closed=:left)
	countPlot = plot(countHist,label=false)
	title!(countPlot,"Total counts displayed = $(sum(countHist.weights))")
	xlabel!(countPlot,"t (ns)")
end

# ╔═╡ 9561be7c-ab38-493a-8b33-cd7fb7aec495
if makePlots && countsPlot
	countEdgesDf   = DataFrame(count_bin_edges = countHist.edges[1])
	countWeightsDf = DataFrame(count_bin_weights = countHist.weights)
	countDataDf    = DataFrame(time_ns = flatγCountTimes)
	
	countEdgesString   = string(collect(CSV.RowWriter(countEdgesDf))...)
	countWeightsString = string(collect(CSV.RowWriter(countWeightsDf))...)
	countDataString    = string(collect(CSV.RowWriter(countDataDf))...)
	
	md"""
	Bin edges:
	$(DownloadButton(countEdgesString,"photon-count_histogram_bin-edges.csv"))
	
	Bin weights:
	$(DownloadButton(countWeightsString,"photon-count_histogram_bin-weights.csv"))
	
	Arrival times:
	$(DownloadButton(countWeightsString,"photon-arrival-time_data.csv"))
	"""
end

# ╔═╡ a62a3b1c-4980-492c-95d5-20681ca343ae
if makePlots && corrPlot
	correlationTimes = map(γCt->singleDeltaTimes(τ,γCt),γCounts)
	flatCorrelationTimes = vcat(correlationTimes...)
end;

# ╔═╡ 25b20ebd-3fc1-4326-9848-1adf5afb768d
if makePlots && corrPlot
	corrHist = fit(Histogram, flatCorrelationTimes,vcat(τ,tres/10)[1+windowPos:windowSize+windowPos],closed=:left)
	plot(corrHist,label=false)
	
end

# ╔═╡ ac647fa5-3294-4328-9b9f-fae54ef54174
if makePlots && corrPlot
	corrEdgesDf   = DataFrame(corr_bin_edges = corrHist.edges[1])
	corrWeightsDf = DataFrame(corr_bin_weights = corrHist.weights)
	corrDataDf    = DataFrame(corr_times = flatCorrelationTimes)
	
	corrEdgesString = string(collect(CSV.RowWriter(corrEdgesDf))...)
	corrWeightsString = string(collect(CSV.RowWriter(corrWeightsDf))...)
	corrDataString = string(collect(CSV.RowWriter(corrDataDf))...)

	md"""
	Bin edges:
	$(DownloadButton(corrEdgesString,"photon-correlation_histogram_bin-edges.csv"))
	
	Bin weights:
	$(DownloadButton(corrWeightsString,"photon-correlation_histogram_bin-weights.csv"))

	Correlation times:
	$(DownloadButton(corrWeightsString,"photon-correlation-tau_data.csv"))
	"""	
	
end

# ╔═╡ 3edfcf46-9f6d-4ff0-8179-f9a7ddb1897f
"""
	function beCount(nbar::Real)

Returns Bose-Einstein distributed counts for average count rate nbar
"""
function beCount(nbar::Real)
	p = 1/(nbar+1)
	fnbar = p*nbar
	f = p*nbar
	s = p
	r = rand()
	count = 0
	while r>s
		count += 1
		p *= f
		s += p
	end
	return count
end

# ╔═╡ 2328359b-8f27-4a18-8d63-0538dbdfcdff
function vectorAvg(someVector::Vector)
	return +(someVector...)/length(someVector)
end

# ╔═╡ 099f95dc-1dbc-48f6-ab3f-0920e86f0c0a
if makePlots && classicalAvgPlots
	efieldtAvg = vectorAvg(efieldt)
	g1τAvg = vectorAvg(g1τ)
	intensitytAvg = vectorAvg(intensityt)
	g2τAvg = vectorAvg(g2τ)
end;

# ╔═╡ 0e50e055-6c57-40d9-96ca-d8c795213999
if makePlots && classicalAvgPlots
	efieldtAvgPlot = plot(times[1+windowPos:windowSize+windowPos],real.(efieldtAvg)[1+windowPos:windowSize+windowPos],label=false,title="Average E(t): N = $(nTrials)")
	xlabel!(efieldtAvgPlot,"t (ns)")

	eCorrAvgPlot = plot(τ[1+windowPos:windowSize+windowPos],g1τAvg[1+windowPos:windowSize+windowPos],label=false,title="Average g1(τ): N = $(nTrials)")
	xlabel!(eCorrAvgPlot,"τ (ns)")

	intensitytAvgPlot = plot(times[1+windowPos:windowSize+windowPos],intensitytAvg[1+windowPos:windowSize+windowPos],label=false,color=2,title="Average I(t): N = $(nTrials)")
	xlabel!(intensitytAvgPlot,"t (ns)")

	iCorrAvgPlot = plot(τ[1+windowPos:windowSize+windowPos],g2τAvg[1+windowPos:windowSize+windowPos],label="data",color=2,title="Average g2(τ): N = $(nTrials)")
	plot!(iCorrAvgPlot,τ[1+windowPos:windowSize+windowPos],transpose(g2τCalc)[1+windowPos:windowSize+windowPos],label="calc",color=3)
	xlabel!(iCorrAvgPlot,"τ (ns)")

	plot(efieldtAvgPlot,intensitytAvgPlot,eCorrAvgPlot,iCorrAvgPlot,layout=(2,2),size=(800,800))

end

# ╔═╡ c023ca90-8bc8-4183-b0b9-9b1a8bfb4b59
if makePlots && classicalAvgPlots
	classicalFieldAvgDf = DataFrame(
		time                 = times,
		realElectricFieldAvg = real.(efieldtAvg),
		imagElectricFieldAvg = imag.(efieldtAvg),
		intensityAvg         = intensitytAvg
		)
	classicalCorrAvgDf = DataFrame(
		tau      = τ,
		g1tauAvg = g1τAvg,
		g2tauAvg = g2τAvg		
		)

	classicalFieldAvgString = string(collect(CSV.RowWriter(classicalFieldAvgDf))...)
	classicalCorrAvgString  = string(collect(CSV.RowWriter(classicalCorrAvgDf))...)
	md"""
	E-field/intensity:
	$(DownloadButton(classicalFieldString,"classical_field_avg.csv"))
	
	Correlations:
	$(DownloadButton(classicalCorrString,"classical_correlations_avg.csv"))
	"""
end

# ╔═╡ 57b0a316-782e-4f85-8a77-8743abc601f1
md"""
## Load Prerequisites
"""

# ╔═╡ Cell order:
# ╟─0a435f6c-9f0a-11eb-114f-fb8c733dafbb
# ╟─5c181bf0-f92e-40ac-9193-62c48927f60f
# ╟─0fe833e8-a359-41e2-ba4a-b8313a222182
# ╟─7ac53f6b-15ee-440d-acae-005dc3acce51
# ╟─c0e353c3-53e7-4419-b8cf-878f6fbc761d
# ╟─402b4a2a-d22f-4ee7-ae3d-ca5fe7fda110
# ╟─89b17676-6496-4d37-929d-3f626e329cdd
# ╟─30f9f516-568b-4e73-8aa2-4f5bb48f9ab4
# ╟─8b607526-47b0-4a8a-aa16-2aee351b590f
# ╟─a3cd4171-5bc5-4988-a5ef-428dd6fa31ab
# ╟─e2fab25a-ccde-44e8-a3bf-e402fc32d2ba
# ╟─57da1ea1-5a5d-49d1-88f2-2624a14f8dff
# ╟─9e013576-7682-4d6c-8caa-e3f2c1034f73
# ╟─b2def348-d4c0-4461-bee6-0373673a8fd8
# ╟─b5342349-5b28-4918-bedb-effa46df10c3
# ╟─1e48204e-23d8-4581-a160-7ba84b83b8d8
# ╟─0e50e055-6c57-40d9-96ca-d8c795213999
# ╟─c023ca90-8bc8-4183-b0b9-9b1a8bfb4b59
# ╟─0f36d9b7-6ba6-402d-9601-53b16e952110
# ╟─92d93e89-1d36-4bf0-b3f8-d4370e08857b
# ╟─9561be7c-ab38-493a-8b33-cd7fb7aec495
# ╟─c54073d5-b58a-46ed-b086-d6347c1dbcda
# ╟─25b20ebd-3fc1-4326-9848-1adf5afb768d
# ╟─ac647fa5-3294-4328-9b9f-fae54ef54174
# ╟─13e70af4-146d-43e7-9fd6-ec8f901185d4
# ╟─807e056a-3c89-46f9-8dc2-5dde04023ce2
# ╠═f2e2b9fd-9d21-410c-89d8-2286b163bb2a
# ╠═21e16aae-1b50-4f31-936a-17525352afbe
# ╟─173ca217-29e6-4c5f-8995-0bd5b62c32dd
# ╠═18f7bc41-751f-4283-8a90-af3ad217b58b
# ╠═39e8d2c6-c8de-4e0a-b531-60dbaadd6aa1
# ╠═099f95dc-1dbc-48f6-ab3f-0920e86f0c0a
# ╟─689bc5a6-2469-41eb-99ce-c300dd73874e
# ╟─f2a04b1f-f61c-45d4-a67a-1609f53d452e
# ╠═045fbba2-5ae7-4a41-9e8e-343a5ab41c7b
# ╟─730e47f5-1c6a-4fd3-a83f-f877b14cde06
# ╠═a62a3b1c-4980-492c-95d5-20681ca343ae
# ╟─282a7d64-8d67-4a10-9aa9-d64be7c5dd54
# ╠═a3e82734-c1d3-4e9b-8603-ce6bc1cb5018
# ╠═d2ec427e-50f8-47d2-9a44-cb0d5c0dcb9b
# ╠═26e08020-6775-431f-b33a-0efaad94189a
# ╠═8175ac47-b261-4942-a70d-f127f22aee74
# ╠═d8c63dd0-d6dd-4234-8515-caf5803af683
# ╠═b237ca1e-1719-420c-820d-67b8e18c7a57
# ╠═3620f708-3eac-4cd8-90db-72ac383b5917
# ╠═48fd09eb-f5a0-4bc4-bdbf-3648366978a1
# ╠═31b92f9f-6be0-48cb-9a7e-bfc5edb79bd2
# ╠═3edfcf46-9f6d-4ff0-8179-f9a7ddb1897f
# ╠═2328359b-8f27-4a18-8d63-0538dbdfcdff
# ╟─57b0a316-782e-4f85-8a77-8743abc601f1
# ╠═f1c9d31f-c565-4e96-8801-22352b014397
# ╠═f6daaa32-fe10-4a20-b30e-ba31b63c6b90
# ╠═338200b2-c10c-4e7d-8a24-bfd32e6757e5
# ╠═525ba450-c635-4bb9-a2f3-aeb12d5c30d3
# ╠═64adb0f4-b695-4be1-9b02-c42562951b39
# ╠═55405eea-1c12-4060-81f5-cbcc082a992c
# ╠═7e392b1a-31a8-4364-ac93-742221501f1c
# ╠═b708f1d2-e8dd-40ee-92fd-cb4ae02a6a9a
# ╠═34a840a1-7050-4719-a81e-c4545551ec7d
# ╠═83920316-034e-4ca1-94fa-5aea4d4d117b
# ╠═8287cc85-37d3-4b53-ab65-0aef5d996a72
# ╠═1196f179-17ff-42d9-a582-1b8498b7c372
# ╠═2d8d5528-d64d-4026-b5a6-93b66e17e609
# ╠═c3929637-0224-4f95-a8ce-b1a159084014
# ╠═2f28edec-8461-4eda-a95f-05f4db5b3a81
