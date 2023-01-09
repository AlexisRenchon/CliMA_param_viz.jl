using WGLMakie, JSServe, SparseArrays
# using GLMakie, SparseArrays
using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
FT = Float64
earth_param_set = create_lsm_parameters(FT)

function figure()
  # Create figure, axis, widgets + layout
  fig = Figure(resolution = (1200, 1200))
  menu = Menu(fig[1,1], options = ["CO2 production", "CO2 diffusivity"])
  ax3D = Axis3(fig[2,1], xlabel = "Soil temperature °C", ylabel = "Soil moisture [m³ m⁻³]")
  axT = Axis(fig[3,1], xlabel = "Soil temperature °C")
  axM = Axis(fig[4,1], xlabel = "Soil moisture [m³ m⁻³]")

  params = SoilCO2ModelParameters{FT}(; earth_param_set = earth_param_set)
  labels = ["$(s)" for s in fieldnames(SoilCO2ModelParameters)[1:end-1]] # without earth_param_set
  ranges = [(val/2 : val/4: val*2, val) for val in [getfield(params, i) for i in fieldnames(SoilCO2ModelParameters)[1:end-1]]]

  sliders = [(label = label, range = range, startvalue = startvalue) for (label, (range, startvalue)) in zip(labels, ranges)]
  sg = SliderGrid(fig[3, 2], sliders..., width = 350, tellheight = false)
  sd = Dict(i => sg.sliders[i].value for i in 1:length(sliders))

  keys = Symbol.(labels) 
  s = collect(values(sort(sd)))
  # ss = @lift($s) # 
  args = (; zip(keys, s)...)
  # args_o = Observable(args) # if we want args to be the Observable

  st = rand(13) # to test with non Observable type
  argss = (; zip(keys, st)...) 
  SoilCO2ModelParameters{FT}(; argss..., earth_param_set = earth_param_set) # this works

  SoilCO2ModelParameters{FT}(; args..., earth_param_set = earth_param_set) # this is equivalent to the second parameters definition below - same error
                                        # how to use @lift($) here?

  @lift(SoilCO2ModelParameters{FT}(; $args..., earth_param_set = earth_param_set)) # MethodError: no method matching iterate(::Observable{Any})

  # TO DO, find a way to not have to write explicitly keyword arguments, so that it can be generalised to any ModelParameters
  parameters = @lift(
                   SoilCO2ModelParameters{FT}(;
                   P_sfc = $(sd[1]),
                   ν = $(sd[2]),
                   θ_a100 = $(sd[3]),
                   D_ref = $(sd[4]),
                   b = $(sd[5]),
                   D_liq = $(sd[6]),
                   α_sx = $(sd[7]),
                   Ea_sx = $(sd[8]),
                   kM_sx = $(sd[9]),
                   kM_o2 = $(sd[10]),
                   O2_a = $(sd[11]),
                   D_oa = $(sd[12]),
                   p_sx = $(sd[13]), 
                   earth_param_set = earth_param_set
                   )
                   )

  parameters = 
                   SoilCO2ModelParameters{FT}(;
                   P_sfc = sd[1],
                   ν = sd[2],
                   θ_a100 = sd[3],
                   D_ref = sd[4],
                   b = sd[5],
                   D_liq = sd[6],
                   α_sx = sd[7],
                   Ea_sx = sd[8],
                   kM_sx = sd[9],
                   kM_o2 = sd[10],
                   O2_a = sd[11],
                   D_oa = sd[12],
                   p_sx = sd[13], 
                   earth_param_set = earth_param_set
                   )
                   



























  fun = [(x, y, p) -> microbe_source(x, y, Csom, p), # need to give a value for Csom
         co2_diffusivity]



#vec = ["d1", "d2", "d3"]
#values = ["foo", "bar", "baz"]
#for i in 1:length(vec)
#    eval(:(vec[i] = "$(values[i])"))
#end


for (i, variable) in enumerate(vec)
    println(i)
    println(variable)
    #eval(:(variable = "$(values[i])"))
end


eval(:(vec[1] = "$(values[1])"))

# Solution below introduces Global variables...
Base.setproperty!(Main, :foo, "bar")

# Solution below could work? Splat a NamedTuple
t = (a=1, b=2, c=3)
function f(x, y, z)
    println("x = $x, y = $y, z = $z")
end
f(t...)

t = [(labels[i] = i) for i = 1:3]

# Somehow use a Dict?

v = ["a", "b", "c"]
t = NamedTuple{v}((i, i) for i in 1:length(v))
t = NamedTuple{v}(Symbol(x) => i for (i, x) in enumerate(v))

# below looks good
keys = (:a, :b, :c)
values = (1, 2, 3)
(; zip(keys, values)...)
  














  ####
  #### Code below needs to be simplified and generalised. 
  ####

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

  x = @lift(mat([10, 40], [0.0, 0.5], 30, fun[$s1], $parameters)[1]) 
  y = @lift(mat([10, 40], [0.0, 0.5], 30, fun[$s1], $parameters)[2])
  z = @lift(mat([10, 40], [0.0, 0.5], 30, fun[$s1], $parameters)[3])

  surface!(ax3D, x, y, z)

  # plot on 2D axis

  x_axT = collect(10:1:40)
  x_axM = collect(0.0:0.0163:0.5)

  y_axT = @lift(fvec(x_axT, $s15, fun[$s1], $parameters))
  y_axM = @lift(fvecM($s14, x_axM, fun[$s1], $parameters))

  lines!(axT, x_axT, y_axT, color = :red, linewidth = 4)
  lines!(axM, x_axM, y_axM, color = :blue, linewidth = 4)

  cM = @lift(repeat([$s15], 31))
  cT = @lift(repeat([$s14], 31))
  lines!(ax3D, x_axT, cM, y_axT, color = :red, linewidth = 4) 
  lines!(ax3D, cT, x_axM, y_axM, color = :blue, linewidth = 4)

  on(sg.sliders[1].value) do val
    autolimits!(ax3D)
    autolimits!(axT)
    autolimits!(axM)
  end

  return fig
end

figure()

