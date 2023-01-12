# Notes 01/12/23
#
# Why do this little web dashboard?
# Help with debugging
# Communication (to funders, between scientists, to public, etc.)
# Web = super accessible (just open an URL), e.g., instant in a meeting if subject randomly comes up
# Web = super accessible (can be integrated in a slidev presentation in 1 line of code)
# etc...
#
# To do: 
# 1. Deal with the porosity issue, again... when soil moisture > porosity, crashes
# 2. Separate the script into functions, write the module
# 3. Register package to Julia registry
# 4. Better layout (can we group multiple block, e.g., SliderGrid and its Label)
# 5. What arguments should be entered by user: AbstractModel, Driver, Driver limit, ... ? 
# 6. Parameters should change dynamically with menu (model) selection
# 7. Little tweaks: better colors, maybe straight lines, output value, datainspector
# 8. Eventually if complete App, 2 menu, 1 to chose AbstractModel, 1 to chose parameterization function
# 9. Currently prototype with SoilCO2Model, need to test if work with others
# 10. Where to host? ClimaCoreMakie? ClimaLSM? New repo?
# 11. Deal with cases where a parameterization function only depends on 1 driver (thus no 3D axis to be done)
# 12. Deal with function that do not use earth_param_set (can use multiple dispatch or if)
# 13. When defining functions, need a fixed driver if more than 2. (e.g., SOM here) This needs to be automatised
# 14. Add units - of parameters and output (soil co2 production)
#
# Questions for Kat:
# 1. Why do we need to do using ClimaLSM.Soil.Biogeochemistry specifically?
# 2. We have 3 drivers explicitely declared with SoilCO2: temperature, moisture, SOM (prescribed for now)
#    Do all AbstractModel declare some drivers? Are they always distinct from parameters? 
#    Are constant parameters always somehow disctinct from variable parameters in their declaration? (e.g., using earth_param_set)
# 3. Can we somehow retrieve the name of parameterization functions of an AbstractModel? My understanding is that we currently have 3 parameterisation functions for SoilCO2, but they don't seem to be contained in anything, so we can't list them somehow
# 4. What automated range of parameters seems to be a good idea? (min:step:max, initial_value), note that it helps to have default values like we did for this (otherwise user would have to provide it) 
# 5. Do you know a better/simpler way to do the mat function below? (compute function on a grid to generate matrix)
# 6. Should drivers have a default value too? (if not, user can enter it)

using WGLMakie, JSServe, SparseArrays
# using GLMakie, SparseArrays
using ClimaLSM # Soil.Biogeochemistry.microbe_source should work without line below
using ClimaLSM.Soil.Biogeochemistry # I am curious as to why this needs to be specifically used even after using ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl")) # this feels a bit weird to me... (why have parameters out of src?)

function figure() # args could be: model::AbstractModel, drivers::?, limits::? (limit of drivers)) 
  # Figure, Axis3 (3D) and Axis (2D), Menu
  fig = Figure(resolution = (1200, 1200))
  labels_d = ["Soil temperature [K]", "Soil moisture [m³ m⁻³]"] # Can be retrieved from PrecribedSoil struct (would be a different struct for prognostic) 
  menu_opt = ["CO2 production", "CO2 diffusivity"] # Should get these name directly from AbtractModel param function names
  menu = Menu(fig[1,1], options = menu_opt); m = menu.selection
  ax3D = Axis3(fig[2,1], xlabel = labels_d[1], ylabel = labels_d[2])
  axT = Axis(fig[3,1], xlabel = labels_d[1])
  axM = Axis(fig[4,1], xlabel = labels_d[2])

  # SliderGrid for parameters
  FT = Float64 
  earth_param_set = create_lsm_parameters(FT)
  params = SoilCO2ModelParameters{FT}(; earth_param_set = earth_param_set)
  labels = ["$(s)" for s in fieldnames(SoilCO2ModelParameters)[1:end-1]] # without earth_param_set
  ranges = [(val/2 : val/4: val*2, val) for val in [getfield(params, i) for i in fieldnames(SoilCO2ModelParameters)[1:end-1]]]
  sliders = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels, ranges)]
  sg = SliderGrid(fig[2, 2], sliders..., width = 350, tellheight = true)
  param_title = Label(fig[1, 2], "Parameters", tellheight = false)

  # SliderGrid for drivers
  startval_d = [293.13, 0.3] # Again, should this values taken from AbstractModel somehow, or given by user with an arg? 
  ranges_d = [(val/2 : val/4: val*2, val) for val in startval_d]
  sliders_d = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels_d, ranges_d)]
  sg_d = SliderGrid(fig[4, 2], sliders_d..., width = 350, tellheight = false)
  param_title_d = Label(fig[3, 2], "Drivers", tellheight = false)

  # Get Observable and their values from SliderGrid
  sd = Dict(i => sg.sliders[i].value for i in 1:length(sliders))
  sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
  sd_d = Dict(i => sg_d.sliders[i].value for i in 1:length(sliders_d))
  sd_v_d = Dict(i => sg_d.sliders[i].value[] for i in 1:length(sliders_d))

  # Create struct of AbtractModel (example SoilCO2ModelParameters for now, but should be an argument of function) from SliderGrid values
  keys = Symbol.(labels) 
  s = collect(values(sort(sd_v)))
  args = (; zip(keys, s)...)
  parameters = Observable(SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set)) # this works

  # Menu select the model to display. x is driver #1 y is driver #2, p is parameters. 
  # This is a tricky part. Some parameterization will depend on 1 or 2 drivers, and sometimes on more stuff (e.g. below, Csom is not a param nor driver)
  # Need to find a way to generalize a solution to all parameterization (could be specified by users in args)
  fun = Dict(menu_opt[1] => (x, y, p) -> microbe_source(x, y, Csom, p), # need to give a value for Csom
             menu_opt[2] => co2_diffusivity) # Notre, name of function should be retrieved from the AbstractModel (how? Kat)

  # Put mat, fvec and fvecM in a separate script
  # These function just compute functions on grids... (in order to plot them) may be a simpler way to do this!
  function mat(x, y, r, fun, params) # x and y, range - 2 values   
    x = collect(range(x[1], length=r, stop=x[2])) # T axis, °C from min to max
    y = collect(range(y[1], length=r, stop=y[2])) # M axis, % from min to max
    X = repeat(1:r, inner=r) # X for DAMM matrix 
    Y = repeat(1:r, outer=r) # Y for DAMM matrix
    X2 = repeat(x, inner=r) # T values to fit DAMM on   
    Y2 = repeat(y, outer=r) # M values to fit DAMM on
    # xy = hcat(X2, Y2) # T and M matrix to create DAMM matrix 
    FMatrix = Matrix(sparse(X, Y, fun.(FT.(X2), FT.(Y2), repeat([params], r*r))))
    return x, y, FMatrix
  end

  function fvec(x, y, fun, params)
    vec = fun.(FT.(x), FT.(repeat([y], 31)), repeat([params], 31))
    return vec
  end

  function fvecM(x, y, fun, params)
    vecM = fun.(FT.(repeat([x], 31)), FT.(y), repeat([params], 31))
    return vecM
  end

  Csom = FT(5.0) # Should be either retrieved from AbstractModel or given by user as an arg
  x = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[1]) 
  y = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[2])
  z = @lift(mat([283, 313], [0.0, 0.5], 30, fun[$m], $parameters)[3])
  surface!(ax3D, x, y, z)

  # plot on 2D axis
  x_axT = collect(283:1:313) # should be ax1 and ax2, drivers can be other than temperature or moisture
  x_axM = collect(0.0:0.0163:0.5)

  # special stuff for driver 1 and driver 2 (e.g., temperature and moisture)
  # so need a separate on(sd[drivers[i]])
  y_axT = @lift(fvec(x_axT, $(sd_d[2]), fun[$m], $parameters)) # fun should be Observable from Menu
  y_axM = @lift(fvecM($(sd_d[1]), x_axM, fun[$m], $parameters))

  lines!(axT, x_axT, y_axT, color = :red, linewidth = 4)
  lines!(axM, x_axM, y_axM, color = :blue, linewidth = 4)

  cM = @lift(repeat([$(sd_d[2])], 31)) # Should output an Observable? See old code from on the on(sd[i]) loop
  cT = @lift(repeat([$(sd_d[1])], 31))
  lines!(ax3D, x_axT, cM, y_axT, color = :red, linewidth = 4) 
  lines!(ax3D, cT, x_axM, y_axM, color = :blue, linewidth = 4)

  for i in 1:length(sd)
    on(sd[i]) do val # update parameters and rescale x and y limits
      # println(sd[9][])
      sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
      s = collect(values(sort(sd_v)))
      args = (; zip(keys, s)...)
      parameters[] = SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set)  # new args
      autolimits!(ax3D)
      autolimits!(axT)
      autolimits!(axM)
    end
  end

  for i in 1:2
    on(sd_d[i]) do val
      autolimits!(axT)
      autolimits!(axM)
    end
  end

  on(menu.selection) do val
    autolimits!(ax3D)
    autolimits!(axT)
    autolimits!(axM)
  end

  return fig
end

figure()

