using WGLMakie, JSServe, SparseArrays
# using GLMakie, SparseArrays
using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
FT = Float64
earth_param_set = create_lsm_parameters(FT)

function figure()
  # Figure, Axis3 (3D) and Axis (2D), Menu
  fig = Figure(resolution = (1200, 1200))
  labels_d = ["Soil temperature °C", "Soil moisture [m³ m⁻³"] # Should this be taken from AbstractModel somehow? or an Arg of this function?
  menu = Menu(fig[1,1], options = ["CO2 production", "CO2 diffusivity"]) # Should get these name directly from AbtractModel param function names
  ax3D = Axis3(fig[2,1], xlabel = labels_d[1], ylabel = labels_d[2])
  axT = Axis(fig[3,1], xlabel = labels_d[1])
  axM = Axis(fig[4,1], xlabel = labels_d[2])

  # SliderGrid for parameters
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
  fun = [(x, y, p) -> microbe_source(x, y, Csom, p), # need to give a value for Csom
         co2_diffusivity] # Notre, name of function should be retrieved from the AbstractModel (how? Kat)

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

  # ###########################
  # Stuff below inside on(sd[i])
  # Should also be generalized to all case of AbstractModel
  # Maybe write a separate script for what is inside on(sd[i])
  # need on(menu.selection) do s (select fun)
  # parameters changes on(sd[i]) ... 
  Csom = FT(5.0) # Should be either retrieved from AbstractModel or given by user as an arg
  x = mat([10, 40], [0.0, 0.5], 30, fun[1], parameters[])[1] # fun should be Observable from Menu
  y = mat([10, 40], [0.0, 0.5], 30, fun[1], parameters[])[2]
  z = mat([10, 40], [0.0, 0.5], 30, fun[1], parameters[])[3]
  surface!(ax3D, x, y, z)

  # plot on 2D axis

  x_axT = collect(10:1:40)
  x_axM = collect(0.0:0.0163:0.5)

  # special stuff for driver 1 and driver 2 (e.g., temperature and moisture)
  # so need a separate on(sd[drivers[i]])
  y_axT = fvec(x_axT, sd_v_d[2], fun[1], parameters[]) # fun should be Observable from Menu
  y_axM = fvecM(sd_v_d[1], x_axM, fun[1], parameters[])

  lines!(axT, x_axT, y_axT, color = :red, linewidth = 4)
  lines!(axM, x_axM, y_axM, color = :blue, linewidth = 4)

  cM = repeat([sd_v_d[2]], 31) # Should output an Observable? See old code from on the on(sd[i]) loop
  cT = repeat([sd_v_d[1]], 31)
  lines!(ax3D, x_axT, cM, y_axT, color = :red, linewidth = 4) 
  lines!(ax3D, cT, x_axM, y_axM, color = :blue, linewidth = 4)

  ########################
  ########################
  on(sd[1]) do val # only slider 1 for now
    println(sd[1][])
    sd_v = Dict(i => sg.sliders[i].value[] for i in 1:length(sliders))
    s = collect(values(sort(sd_v)))
    args = (; zip(keys, s)...)
    parameters[] = SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set)  # new args
    autolimits!(ax3D)
    autolimits!(axT)
    autolimits!(axM)
  end

  return fig
end

figure()

