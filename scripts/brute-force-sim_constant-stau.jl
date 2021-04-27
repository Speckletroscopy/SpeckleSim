### A Pluto.jl notebook ###
# v0.14.3

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

# ╔═╡ 42c0a736-b762-40e3-a4fb-46b6303f0cff
Pkg.add("FFTW")

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

# ╔═╡ 258b3c88-b1ca-49b2-b414-8ce8db178724
using FFTW

# ╔═╡ 0a435f6c-9f0a-11eb-114f-fb8c733dafbb
md"# Brute Force Simulation"

# ╔═╡ 5c181bf0-f92e-40ac-9193-62c48927f60f
md"""
## Simulation Configuration
"""

# ╔═╡ 0fe833e8-a359-41e2-ba4a-b8313a222182
md"""
### Classical field parameters

$(@bind bigN NumberField(1:10000;default=100))

Number of emitting atoms


$(@bind temp NumberField(293:3000;default=400))

Temperature (K)

$(@bind tres NumberField(10:500))

Time resolution (ps)

$(@bind tmax NumberField(1.0:0.1:50.0))

Maximum time (ns)

"""

# ╔═╡ a8b59a40-e281-4e59-97ae-c523be21e9eb
md"""
$(@bind shift Slider(1:30,show_value=true ))

Line separation (GHz)
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
#### Run Simulation $\rightarrow$ [$(@bind sim CheckBox())]
"""

# ╔═╡ 659c5d8a-6722-4c74-834a-f79bc4312e8e
if sim
	md"""
	Make instances: $(@bind makeInstances CheckBox())
	"""
end

# ╔═╡ 30f9f516-568b-4e73-8aa2-4f5bb48f9ab4
md"""
## Plots
"""

# ╔═╡ 2f0b2664-8969-4712-b26d-fbaa06b36ae3
if sim
md"""
### Plots from simulated photon counts
"""
end

# ╔═╡ d512f123-088b-4a65-965e-af512147249d
if sim
	md"""
	### Plots using classical field
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

# ╔═╡ 93a1cad0-6b86-40a0-b122-19b854f51174
# define beamsplitter
begin
	bs = [1,1]
	bs = bs/norm(bs)
end;

# ╔═╡ 282a7d64-8d67-4a10-9aa9-d64be7c5dd54
md"""
## Functions
"""

# ╔═╡ c84cf75f-17da-4d54-9e97-9bf69968732c
md"""
### Electric field parameters and generation
"""

# ╔═╡ 26e08020-6775-431f-b33a-0efaad94189a
begin
	"""
		eFieldParams(Em::Vector,ωm::Vector,ω0::Number,σ::Number,rng::MersenneTwister)

	Static parameters for the electric field
	"""
	struct eFieldParams
		Em::Vector # field magnitudes
		ωm::Vector # central frequencies of lines
		ω0::Number # central frequency for composite line
		σ::Number # Doppler broadening
		rng::MersenneTwister
		function eFieldParams(Em::Vector,ωm::Vector,ω0::Number,σ::Number,rng::MersenneTwister)
			@assert length(Em) == length(ωm) "Length of field magnitude and frequency vectors must match"
			new(
				convert(Vector{Complex},Em),
				convert(Vector{Real},ωm),
				convert(Real,ω0),
				convert(Real,σ),
				rng
			)
		end
	end

	
	"""
		eFieldParams(Em::Vector,ωm::Vector,σ::Number,seed::Integer=-1)

	Static parameters for the electric field
	"""
	function eFieldParams(Em::Vector,ωm::Vector,σ::Number,seed::Integer = -1)
		if seed == -1
			rng = MersenneTwister()
		else
			rng = MersenneTwister(seed)
		end
		return eFieldParams(Em,ωm,ωm[1],σ,rng)
	end
end

# ╔═╡ e51b2eac-2de8-4f30-bc11-26c320bd6afb
"""
	function seed!(params::eFieldParams,seed::Integer)

Reset the seed for the random number generator in params
"""
function seed!(params::eFieldParams,seed::Integer)
	Random.seed!(params.rng,seed)
end

# ╔═╡ 904ec626-6678-4f19-8247-e404fe9d278b
begin
	"""
	Container holding frequencies and phases for one realization of the electric field.
	"""
	struct eFieldInstance
		ωn::Vector
		ϕmn::Matrix
		params::eFieldParams

		function eFieldInstance(ωn::Vector,ϕmn::Matrix,params::eFieldParams)
			@assert (size(ϕmn)[1] == length(params.ωm) && size(ϕmn)[2] == length(ωn)) "Phase array must have shape (m,n)"
			new(ωn,ϕmn,params)
		end
	end

	"""
		eFieldInstance(n::Integer,params::eFieldParams)
	
	Generate frequencies and phases for the electric field from n emitting atoms and return an eFieldInstance object.
	"""
	function eFieldInstance(n::Integer,params::eFieldParams)
		ωDist = Normal(params.ω0,params.σ)
		ωn    = rand(params.rng,ωDist,n)
		ϕmn   = 2*π*rand(params.rng,Float64,(length(params.ωm),n))
		return eFieldInstance(ωn,ϕmn,params)
	end
end

# ╔═╡ 8175ac47-b261-4942-a70d-f127f22aee74
"""
	function electricField(t::Real,params::eFieldParams)

Returns the electric field value at time t
"""
function electricField(t::Real,instance::eFieldInstance)
	Δm = instance.params.ωm .- instance.params.ω0
	# generate frequencies
	ωmn = transpose(instance.ωn) .+ Δm
	# add the phase
	exponentmn = -im*(t*ωmn+instance.ϕmn)
	# put them in the exponent
	enm = exp.(exponentmn)
	# multiply by the field magnitude
	fieldnm = instance.params.Em .* enm
	# add it all together
	return sum(ivec(fieldnm))
end

# ╔═╡ 58b00c01-9dca-4220-97a6-5519d7672bd5
"""
	function σTemp(ω0::Real,temp::Number)

Returns the standard deviation due to Doppler broadening of frequency ω0 at temp
"""
function σTemp(ω0::Real,temp::Number)
	kbOverMhC2 = 9.178e-14;
	return sqrt(kbOverMhC2*temp)*ω0
end

# ╔═╡ f2e2b9fd-9d21-410c-89d8-2286b163bb2a
# these parameters are needed for all simulations
begin
	# Balmer-α lines specified here
	# ωM = 2*π*[456811.0, 456812.0]
	ωM = 2*π*[456811.0, 456811.0+shift]
	# magnitude of each line
	mag = convert(Vector{ComplexF64},ones(length(ωM)))
	# Doppler broadening
	σDopp = σTemp(ωM[1],temp)
	# Store parameters
	eParams = eFieldParams(mag,ωM,σDopp)
	
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
if sim
	md"""
	Size: $(@bind windowSize Slider(1:window,default=window))
	"""
end

# ╔═╡ e2fab25a-ccde-44e8-a3bf-e402fc32d2ba
if sim
	md"""
	Position: $(@bind windowPos Slider(0:(window-windowSize)))
	"""
end

# ╔═╡ b61954a4-a0ba-4f88-90af-0baa98d1325c
if sim
	windowRange = 1+windowPos:windowSize+windowPos
end;

# ╔═╡ 8b607526-47b0-4a8a-aa16-2aee351b590f
if sim
	md"""
	**Window:**

	Window size (ns) = $(tres*windowSize/1000)

	Window position (ns) = $(windowPos*tres/1000)
	"""
end

# ╔═╡ 2f81c1e8-b9bd-4b41-8d9f-f54f964518cf
md"""
### Classical Intensity
"""

# ╔═╡ 5b656847-12d5-442d-8c2a-bbdbe3de0ab7
begin
	function intensity(eFieldT::Number)
		return real(eFieldT*conj(eFieldT))
	end

	function intensity(t::Real,instance::eFieldInstance)
		eFieldT = electricField(t,instance)
		return intensity(eFieldT)
	end
end

# ╔═╡ b409aae3-7ea5-4278-b56a-463718b40426
md"""
### Classical Correlations
"""

# ╔═╡ a8a6e7d0-d871-40e4-ac37-9685cd9e001b
md"""
- Our calculated version

$g^{(2)}(\tau) =1-\frac{1}{N}\frac{\sum_{m=1}^M|\mathcal{E}|_m^4}{\left(\sum_{m=1}^M|\mathcal{E}_m|^2\right)^2}+\left|\frac{\sum_{m=1}^M|\mathcal{E}_m|^2e^{-i\Delta_m\tau}}{\sum_{m=1}^M|\mathcal{E}|_m^2}\right|^2\frac{\langle S(\tau)\rangle_\omega}{N^2}$
"""

# ╔═╡ 4bf8740c-b840-4e3f-9413-27490cc9f79b
begin
	function stauAvg(τ::Number,params::eFieldParams,n::Integer)
		return stauAvg(τ,params.σ,n)
	end
	
	function stauAvg(τ::Number,σ::Number,n::Integer)
		term1 = n
		term2 = n*(n-1)
		term2 *= exp(-σ^2*τ^2)
		return term1 + term2
	end
end

# ╔═╡ 9df57654-c599-4313-87be-3e1beeb82448
begin
	function stauVar(τ::Number,params::eFieldParams,n::Integer)
		return stauVar(τ,params.σ,n)
	end
	
	function stauVar(τ::Number,σ::Number,n::Integer)
		στ2 = σ^2*τ^2

		prod1 = 8*n*(n-1)
		prod1 *= exp(-2*στ2)

		prod2 = n-1+cosh(στ2)

		prod3 = sinh(στ2/2)^2

		return prod1*prod2*prod3
	end
end

# ╔═╡ fd228d30-48e2-4b8f-96a4-1c846743686a
begin
	function stau(τ::Number,instance::eFieldInstance)
		return stau(τ,instance.ωn)
	end
	
	function stau(τ::Number,ωn::Vector)
		terms = exp.(-im*τ*ωn)
		sumterms = sum(terms)
		return real(sumterms*conj(sumterms))
	end
end

# ╔═╡ 637c1989-e452-4eea-9a11-3d94f3e2f6f5
function g2Calc(τ::Number,n::Integer,params::eFieldParams)
	em2 = real.(params.Em .* conj.(params.Em))
	em4 = em2 .* em2
	sumEm2 = sum(em2)
	sumEm4 = sum(em4)
	
	g2τ = 1.0
	
	term2 = -sumEm4/(n*sumEm2^2)
	
	g2τ += term2
	
	Δm = params.ωm .- params.ω0
	term3 = sum(em2 .* exp.(-im*τ*Δm))/sumEm2
	term3 *= conj(term3)
	term3 = real(term3)
	term3 *= stauAvg(τ,params,n)/n^2
	
	return g2τ+term3
end

# ╔═╡ 17d19e10-604f-4a5f-9ed9-73bca2f68d4b
md"""
### Fourier transform related functions for classical case
"""

# ╔═╡ 664ede52-c0c5-4591-b94b-e11abd2de1ff
function fftPositiveFreq(ft,f)
	@assert length(ft) == length(f) "Fourier transform and frequency vectors much have matching length"
	fftRange = f .>= 0
	return (ft[fftRange],f[fftRange])
end

# ╔═╡ 3c7feb92-1223-4428-a86c-0b46960df59f
function g2ClassicalFFT(g2τ::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}
	@assert length(τ) == length(g2τ) "Time and coherence arrays must have the same length"
	g2τPrep = g2τ[1+cuts[1]:end-cuts[2]]
	g2τPrepMean = mean(g2τPrep)
	g2τPrep = g2τPrep .- g2τPrepMean
	g2τFFT = fft(g2τPrep)
end

# ╔═╡ fe3043d0-5417-4a0a-90bc-03f82186edd9
function g2ClassicalFreq(tres::Real,τ::Vector,cuts::Tuple{T,T}) where {T<:Integer}
	τlen = length(τ) - sum(cuts)
	freqFFT = fftfreq(τlen,1/tres)
end

# ╔═╡ 9c2a61f6-dc60-4d35-bfaf-45aad4de2905
md"""
### General correlation functions
"""

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

# ╔═╡ 361cec13-e752-4166-8cd3-d74cb532dd56
md"""
### Photon counts from intensity
"""

# ╔═╡ 48fd09eb-f5a0-4bc4-bdbf-3648366978a1
"""
	function γIntensity(intensity::Vector,nbar::Real)

Calculates the photon count rate in each bin of an intensity histogram
"""
function γIntensity(intensity::Vector,nbar::Real)
	nintensity = intensity/sum(intensity)
	return nbar*nintensity
end

# ╔═╡ 933ae343-2226-4ce0-a9e7-9bc581b07de2
md"""
### Photon correlations
"""

# ╔═╡ c2b7737e-faab-4ed4-9b03-0e1286c635cd
"""
	function autocorrTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function autocorrTimes(τ::Vector,γCounts::Vector)
	return τ[map(i->autocorrelate(γCounts,i,length(τ)) > 0 ? true : false,collect(0:length(τ)-1))]
end

# ╔═╡ 43c5fdbe-823d-4b1e-8801-e41e467ffeef
"""
	function corrTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function corrTimes(τ::Vector,γCounts1::Vector,γCounts2::Vector)
	@assert length(γCounts1) == length(γCounts2) "Count vectors must have the same length"
	return τ[map(i->correlate(γCounts1,γCounts2,i,length(τ)) > 0 ? true : false,collect(0:length(τ)-1))]
end

# ╔═╡ 0abe72e2-92c4-4752-856f-b5975ab4c5e3
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

# ╔═╡ 727d4c4c-d374-4b3b-8f4b-8e73b27edc38
md"""
### Count generators
"""

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

# ╔═╡ 18f7bc41-751f-4283-8a90-af3ad217b58b
# This is the main calculation cell for the simulation
if sim	
	a = 1
	# low cut on τ
	τstart =100
	# high cut on τ
	τend = 0

	# generate frequency bins of fourier transform
	allfreqs = g2ClassicalFreq(tres*1e-3,τ,(τstart,τend))

	# array to accumulate g2(τ) Fourier transform
	g2τ1FFTsum = zeros(ComplexF64,length(τ)-τstart)
	
	# array to accumulate all photon counts 
	totCounts = zeros(Integer,length(times))
	
	# concatenation of all photon correlation times
	allCorrTimes = Vector{Float64}[]
	
	### Calculate initial instance of electric field ###
	
	# instantiate electric field 
	eInstance1 = eFieldInstance(bigN,eParams)

	# calculate electric field vs time
	eFieldT1 = map(t->electricField(t,eInstance1),times)

	# apply beam splitter
	eFieldT1Beam = bs[1] * eFieldT1

	# calculate intensity
	intensity1Beam = intensity.(eFieldT1Beam)

	# calculate classical g2τ
	g2τ1Norm = mean(intensity1Beam)^2
	g2τ1 = map(i->autocorrelate(intensity1Beam,i,window),collect(0:window-1))/g2τ1Norm

	# calculate Fourier transform of this instance
	g2τ1FFT = g2ClassicalFFT(g2τ1,(τstart,τend))		

	# accumulate Fourier transform sum
	g2τ1FFTsum .+= g2τ1FFT		

	# number of instances
	trials = convert(Integer,ceil(ntot/nbar))
	
	# loop over instances
	for n = 1:trials
		
		# calculate the average photon counts in each time bin from the beam intensity and the overall average photon count rate
		γavg1Beam = γIntensity(intensity1Beam,nbar/2)
		
		# generate counts for each beam
		γcounts1Beam1 = poissonCount.(γavg1Beam)
		γcounts1Beam2 = poissonCount.(γavg1Beam)
		
		# concatenate the correlations from this run into all others
		allCorrTimes = vcat(allCorrTimes,corrTimes(τ,γcounts1Beam1,γcounts1Beam2))
		
		# add counts from both beams
		totCounts .+= γcounts1Beam1
		totCounts .+= γcounts1Beam2
		
		# reinstantiate if desired
		if makeInstances
			# reinstantiate electric field 
			eInstance1 = eFieldInstance(bigN,eParams)

			# calculate electric field vs time
			eFieldT1 = map(t->electricField(t,eInstance1),times)

			# apply beam splitter
			eFieldT1Beam = bs[1] * eFieldT1

			# calculate intensity
			intensity1Beam = intensity.(eFieldT1Beam)

			# calculate classical g2τ
			g2τ1Norm = mean(intensity1Beam)^2
			g2τ1 = map(i->autocorrelate(intensity1Beam,i,window),collect(0:window-1))/g2τ1Norm

			# calculate Fourier transform of this instance
			g2τ1FFT = g2ClassicalFFT(g2τ1,(τstart,τend))		

			# accumulate Fourier transform sum
			g2τ1FFTsum .+= g2τ1FFT		
		end
	end
	
	g2τ1FFTsinglePos,allfreqsPos = fftPositiveFreq(g2τ1FFT,allfreqs)
	g2τ1FFTsumPos,_ = fftPositiveFreq(g2τ1FFTsum,allfreqs)
	
end;

# ╔═╡ 4a6763e9-46c3-470e-8712-2e4aa65493cf
if sim
	γCountsPlot = bar(times[windowRange],totCounts[windowRange],label=false)
	title!(γCountsPlot,"Total photon counts vs time")
	
end

# ╔═╡ 5fdf6412-ab5e-434b-b5df-3c16fc98216d
# make electric field plot
if sim
	eFieldPlot = plot(times[windowRange],real.(eFieldT1Beam[windowRange]),label=false)
	xlabel!(eFieldPlot,"t (ns)")
	title!(eFieldPlot,"E(t)")
end

# ╔═╡ dde1cc54-1ae5-4c0f-8762-b1cf3e2f30e9
# make intensity plot
if sim
	intensityPlot = plot(times[windowRange],intensity1Beam[windowRange],label=false)
	xlabel!(intensityPlot,"t (ns)")
	title!(intensityPlot,"I(t)")
end

# ╔═╡ 0591df60-f054-4136-819b-909bad833ed6
if sim
	g2τPlot = plot(τ[windowRange],g2τ1[windowRange],label=false)
	xlabel!(g2τPlot,"τ (ns)")
	title!(g2τPlot,"g2(τ)")
end

# ╔═╡ 94f06e2a-bf0f-4d78-b3a9-83e7865d05a1
if sim
	g2τFFTplot = plot(allfreqsPos,abs.(g2τ1FFTsinglePos),label=false)
end

# ╔═╡ a2d886d7-f030-4d38-860b-399e501b24fc
if sim && makeInstances
	g2τFFTsumPlot = plot(allfreqsPos,abs.(g2τ1FFTsumPos),label=false)
end

# ╔═╡ dd5ba54e-5049-414c-af9f-f8cb1925b662
# make photon count correlation histogram and take fourier transform
if sim
	
	# bin correlation times into a histogram
	corrHist = fit(Histogram,allCorrTimes,vcat(τ,τ[end]+tres*1e-3),closed=:left)
	
	# normalize the histogram
	normCorr = corrHist.weights/sum(corrHist.weights)
	
	# subtract mean from histogram to reduce constant term
	normAvgCorr = normCorr .- mean(normCorr)
	
	# take the fourier transform
	γfft = fft(normAvgCorr)
	
	# recover fourier transform frequencies
	γfreqs = fftfreq(length(γfft),1/(tres*1e-3))
	
	# select positive frequencies for plotting
	γfftPos,γfreqsPos = fftPositiveFreq(γfft,γfreqs)
end;

# ╔═╡ 34a9830b-8b4a-441d-a1db-20ae72565ef5
if sim
	md"""
	Size: $(@bind fwindowSize Slider(1:length(γfreqsPos),default=length(γfreqsPos)-10))
	"""
end

# ╔═╡ 92980437-333e-4205-9b53-9de074a849d4
if sim
	md"""
	Position: $(@bind fwindowPos Slider(0:(length(γfreqsPos)-fwindowSize),default=5))
	"""
end

# ╔═╡ 730295a5-5217-4a80-a327-b0d61a88165a
if sim
	fwindow = 1+fwindowPos:fwindowPos+fwindowSize
end;

# ╔═╡ 9f1f5473-1acd-4d13-a7aa-e88dcdf2420a
if sim
	γCorrFreqPlot = plot(γfreqsPos[fwindow],abs.(γfftPos)[fwindow],label = false)
	xlabel!(γCorrFreqPlot,"frequency (GHz)")
	title!(γCorrFreqPlot,"Fourier transform of photon correlations")
end

# ╔═╡ 50c3ba83-89fd-4211-a3d1-1439e6a7cb5e
if sim
	γCorrPlot = plot(corrHist,label=false)
	xlabel!(γCorrPlot,"τ (ns)")
	title!(γCorrPlot,"Correlated counts vs correlation time")
	
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

# ╔═╡ 72cbac0a-ae8a-4296-ba78-ab5696539a5a
md"""
### Helping functions
"""

# ╔═╡ 2328359b-8f27-4a18-8d63-0538dbdfcdff
function vectorAvg(someVector::Vector)
	return +(someVector...)/length(someVector)
end

# ╔═╡ 57b0a316-782e-4f85-8a77-8743abc601f1
md"""
## Load Prerequisites
"""

# ╔═╡ Cell order:
# ╟─0a435f6c-9f0a-11eb-114f-fb8c733dafbb
# ╟─5c181bf0-f92e-40ac-9193-62c48927f60f
# ╟─0fe833e8-a359-41e2-ba4a-b8313a222182
# ╟─a8b59a40-e281-4e59-97ae-c523be21e9eb
# ╟─7ac53f6b-15ee-440d-acae-005dc3acce51
# ╟─c0e353c3-53e7-4419-b8cf-878f6fbc761d
# ╟─659c5d8a-6722-4c74-834a-f79bc4312e8e
# ╟─402b4a2a-d22f-4ee7-ae3d-ca5fe7fda110
# ╟─89b17676-6496-4d37-929d-3f626e329cdd
# ╟─30f9f516-568b-4e73-8aa2-4f5bb48f9ab4
# ╟─2f0b2664-8969-4712-b26d-fbaa06b36ae3
# ╟─730295a5-5217-4a80-a327-b0d61a88165a
# ╟─34a9830b-8b4a-441d-a1db-20ae72565ef5
# ╟─92980437-333e-4205-9b53-9de074a849d4
# ╠═9f1f5473-1acd-4d13-a7aa-e88dcdf2420a
# ╠═50c3ba83-89fd-4211-a3d1-1439e6a7cb5e
# ╠═4a6763e9-46c3-470e-8712-2e4aa65493cf
# ╟─d512f123-088b-4a65-965e-af512147249d
# ╟─b61954a4-a0ba-4f88-90af-0baa98d1325c
# ╟─8b607526-47b0-4a8a-aa16-2aee351b590f
# ╟─a3cd4171-5bc5-4988-a5ef-428dd6fa31ab
# ╟─e2fab25a-ccde-44e8-a3bf-e402fc32d2ba
# ╟─5fdf6412-ab5e-434b-b5df-3c16fc98216d
# ╟─dde1cc54-1ae5-4c0f-8762-b1cf3e2f30e9
# ╠═0591df60-f054-4136-819b-909bad833ed6
# ╠═94f06e2a-bf0f-4d78-b3a9-83e7865d05a1
# ╠═a2d886d7-f030-4d38-860b-399e501b24fc
# ╟─13e70af4-146d-43e7-9fd6-ec8f901185d4
# ╟─807e056a-3c89-46f9-8dc2-5dde04023ce2
# ╠═f2e2b9fd-9d21-410c-89d8-2286b163bb2a
# ╠═93a1cad0-6b86-40a0-b122-19b854f51174
# ╠═18f7bc41-751f-4283-8a90-af3ad217b58b
# ╠═dd5ba54e-5049-414c-af9f-f8cb1925b662
# ╟─282a7d64-8d67-4a10-9aa9-d64be7c5dd54
# ╟─c84cf75f-17da-4d54-9e97-9bf69968732c
# ╠═26e08020-6775-431f-b33a-0efaad94189a
# ╠═e51b2eac-2de8-4f30-bc11-26c320bd6afb
# ╠═904ec626-6678-4f19-8247-e404fe9d278b
# ╠═8175ac47-b261-4942-a70d-f127f22aee74
# ╠═58b00c01-9dca-4220-97a6-5519d7672bd5
# ╟─2f81c1e8-b9bd-4b41-8d9f-f54f964518cf
# ╠═5b656847-12d5-442d-8c2a-bbdbe3de0ab7
# ╟─b409aae3-7ea5-4278-b56a-463718b40426
# ╟─a8a6e7d0-d871-40e4-ac37-9685cd9e001b
# ╠═4bf8740c-b840-4e3f-9413-27490cc9f79b
# ╠═9df57654-c599-4313-87be-3e1beeb82448
# ╠═fd228d30-48e2-4b8f-96a4-1c846743686a
# ╠═637c1989-e452-4eea-9a11-3d94f3e2f6f5
# ╟─17d19e10-604f-4a5f-9ed9-73bca2f68d4b
# ╠═664ede52-c0c5-4591-b94b-e11abd2de1ff
# ╠═3c7feb92-1223-4428-a86c-0b46960df59f
# ╠═fe3043d0-5417-4a0a-90bc-03f82186edd9
# ╟─9c2a61f6-dc60-4d35-bfaf-45aad4de2905
# ╠═b237ca1e-1719-420c-820d-67b8e18c7a57
# ╠═3620f708-3eac-4cd8-90db-72ac383b5917
# ╟─361cec13-e752-4166-8cd3-d74cb532dd56
# ╠═48fd09eb-f5a0-4bc4-bdbf-3648366978a1
# ╟─933ae343-2226-4ce0-a9e7-9bc581b07de2
# ╠═c2b7737e-faab-4ed4-9b03-0e1286c635cd
# ╠═43c5fdbe-823d-4b1e-8801-e41e467ffeef
# ╠═0abe72e2-92c4-4752-856f-b5975ab4c5e3
# ╟─727d4c4c-d374-4b3b-8f4b-8e73b27edc38
# ╠═31b92f9f-6be0-48cb-9a7e-bfc5edb79bd2
# ╠═3edfcf46-9f6d-4ff0-8179-f9a7ddb1897f
# ╟─72cbac0a-ae8a-4296-ba78-ab5696539a5a
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
# ╠═42c0a736-b762-40e3-a4fb-46b6303f0cff
# ╠═34a840a1-7050-4719-a81e-c4545551ec7d
# ╠═83920316-034e-4ca1-94fa-5aea4d4d117b
# ╠═8287cc85-37d3-4b53-ab65-0aef5d996a72
# ╠═1196f179-17ff-42d9-a582-1b8498b7c372
# ╠═2d8d5528-d64d-4026-b5a6-93b66e17e609
# ╠═c3929637-0224-4f95-a8ce-b1a159084014
# ╠═2f28edec-8461-4eda-a95f-05f4db5b3a81
# ╠═258b3c88-b1ca-49b2-b414-8ce8db178724
