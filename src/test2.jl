using WGLMakie, JSServe, SparseArrays

using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry

# Why is this not is ClimaLSM src?
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
FT = Float64
earth_param_set = create_lsm_parameters(FT)

# Makie
function figure()
  # Create figure, axis, widgets + layout
  fig = Figure(resolution = (1200, 1200))
  menu = Menu(fig[1,1], options = ["CO2 production", "CO2 diffusivity"])
  ax3D = Axis3(fig[2,1], xlabel = "Soil temperature °C", ylabel = "Soil moisture [m³ m⁻³]")
  axT = Axis(fig[3,1], xlabel = "Soil temperature °C")
  axM = Axis(fig[4,1], xlabel = "Soil moisture [m³ m⁻³]")

  params = SoilCO2ModelParameters{FT}(; earth_param_set = earth_param_set)
  labels = ["$(s)" for s in fieldnames(SoilCO2ModelParameters)[1:end-1]] # without earth_param_set
  ranges = [[((getfield(params, i))/2 : (getfield(params, i))/4: (getfield(params, i))*2, getfield(params, i))] for i in fieldnames(SoilCO2ModelParameters)[1:end-1]]

  sliders = []
  for (label, (range, startvalue)) in zip(labels, ranges)
    push!(sliders, (label = label, range = range, startvalue = startvalue))
  end

  sg = SliderGrid(fig[2, 2], sliders..., width = 350, tellheight = false)

  s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15 = sg.sliders[:] # need to test

  return fig
end

figure()












