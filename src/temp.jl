using ClimaLSM
using ClimaLSM.Soil.Biogeochemistry

# works with either GLMakie or WGLMakie

using WGLMakie, SparseArrays

function testing()
fig = Figure(resolution = (1200, 1200))
ax3D = Axis3(fig[1,1], xlabel = "Soil temperature °C", ylabel = "Soil moisture [m³ m⁻³]")
axT = Axis(fig[2,1], xlabel = "Soil temperature °C")
axM = Axis(fig[3,1], xlabel = "Soil moisture [m³ m⁻³]")

FT = Float32

# Parameters should be supplied in m/kg/s (Pa... etc)
P_sfc = FT(101e3) # 1 [Pa] pressure just above the soil surface at time t
Rb = FT(6e-5) # 5 [kg m-3] total root biomass C in soil column
α1r = FT(11.65) # 6 [-]
α2r = FT(20.7) # 7 [-]
α3r = FT(-164.2) # 8 [-]
Vb = FT(0.0015) # 11 [kg C m-3 s-1] value of Vmax at 10 °C and mean environmental conditions
α1m = FT(14.05) # 12 [-]
α2m = FT(11.05) # 13 [-]
α3m = FT(-87.6) # 14 [-]
Km = FT(10e-5) # 15 [kg C m-3 s-1] Michaelis-Menten half saturation constant
CUE = FT(0.8) # 16 [kg C kg-1 C-1] microbial carbon use efficiency
soluble_fraction = FT(0.004) # 17 [-] fraction of soil organic C that is soluble
D_liq = FT(3.17) # 18 [-] Diffusivity of soil C substrate in liquid
Estar = FT(324.6) # 19 [Kelvin] temperature sensitivity parameter
T_ref = FT(273.15) # 20 [Kelvin] temperature sensitivity-related parameter
α4 = FT(4.7) # 21 [-]
T_ref_soil = FT(283.15) # 22 [Kelvin] ref temperature for other param e.g., Rb
α5 = FT(4.547) # 23 [-]
ν = FT(0.556) # 26 [m3 m-3] soil porosity
θ_a100 = FT(0.1846) # 25 air filled porosity at soil water potential of -100 cm H2O
D_ref = FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
P_ref = FT(101325) # 28 [Pa] standard pressure
b = FT(4.547) # 29 [-] parameter related to pore size distribution

# Prognostic variables
T_soil = FT(20)
θ_l = FT(0.3)
θ_i = FT(0.0)
θ_w = θ_l + θ_i
θ_ant_roots = FT(0.3)
θ_ant_microbe = FT(0.3)
T_ant_soil = FT(303)
Cr = FT(10.0) # 2 [kg C m-3] Cr is root biomass carbon, see figure S5
Csom = FT(5.0) # 3 [kg C m-3] soil organic C content at depth z
Cmic = FT(1.0) # 4 [kg C m-3] Microbial C pool, ~ 1 % of Csom at DETECT site

sg = SliderGrid( # probably not implemented in WGLMakie, but convenient in GLMakie
    fig[2, 2],
    (label = "Function", range = 1:1:3, startvalue = 1),
    (label = "Rb", range = Rb - Rb*2 : Rb/2 : Rb + Rb*2, startvalue = Rb),
    (label = "Vb", range = Vb - Vb*2 : Vb/2 : Vb + Vb*2, startvalue = Vb),
    (label = "Km", range = Km - Km*2 : Km/2 : Km + Km*2, startvalue = Km),
    (label = "CUE", range = CUE - CUE*2 : CUE/2 : CUE + CUE*2, startvalue = CUE),
    (label = "Estar", range = Estar - Estar*2 : Estar/2 : Estar + Estar*2, startvalue = Estar),
    (label = "ν", range = ν - ν*2 : ν/2 : ν + ν*2, startvalue = ν),
    (label = "T_ant_soil", range = T_ant_soil - T_ant_soil*2 : T_ant_soil/2 : T_ant_soil + T_ant_soil*2, startvalue = T_ant_soil),
    (label = "θ_ant_roots", range = θ_ant_roots - θ_ant_roots*2 : θ_ant_roots/2 : θ_ant_roots + θ_ant_roots*2, startvalue = θ_ant_roots),
    (label = "θ_ant_microbe", range = θ_ant_microbe - θ_ant_microbe*2 : θ_ant_microbe/2 : θ_ant_microbe + θ_ant_microbe*2, startvalue = θ_ant_microbe),
    (label = "Cr", range = Cr - Cr*2 : Cr/2 : Cr + Cr*2, startvalue = Cr),
    (label = "Csom", range = Csom - Csom*2 : Csom/2 : Csom + Csom*2, startvalue = Csom),
    (label = "Cmic", range = Cmic - Cmic*2 : Cmic/2 : Cmic + Cmic*2, startvalue = Cmic),
    (label = "T_soil", range = 10:1:40, startvalue = T_soil),
    (label = "θ_l", range = 0.0:0.0163:0.5, startvalue = θ_l),
    width = 350,
    tellheight = false)

s1 = sg.sliders[1].value # fun
s2 = sg.sliders[2].value
s3 = sg.sliders[3].value
s4 = sg.sliders[4].value
s5 = sg.sliders[5].value
s6 = sg.sliders[6].value
s7 = sg.sliders[7].value
s8 = sg.sliders[8].value # T_ant_soil
s9 = sg.sliders[9].value
s10 = sg.sliders[10].value 
s11 = sg.sliders[11].value
s12 = sg.sliders[12].value
s13 = sg.sliders[13].value
s14 = sg.sliders[14].value
s15 = sg.sliders[15].value

parameters = @lift(
                   SoilCO2ModelParameters(;
                   P_sfc = P_sfc,
                   Rb = $s2,
                   α1r = α1r,
                   α2r = α2r,
                   α3r = α3r,
                   Vb = $s3,
                   α1m = α1m,
                   α2m = α2m,
                   α3m = α3m,
                   Km = $s4,
                   CUE = $s5,
                   soluble_fraction = soluble_fraction,
                   D_liq = D_liq,
                   Estar = $s6,
                   T_ref = T_ref,
                   α4 = α4,
                   T_ref_soil = T_ref_soil,
                   α5 = α5,
                   ν = $s7, 
                   θ_a100 = θ_a100,
                   D_ref = D_ref,
                   P_ref = P_ref,
                   b = b,
                   )
                   )

fun = [
       #volumetric_air_content, # θ_a = f(θ_w; ν) 
       #co2_diffusivity, # D = f(T_soil, θ_w; θ_a, P_sfc, D_ref, T_ref, P_ref, θ_a100, b, ν)
       #root_source_moisture_coeff, # coeff = f(θ_l, θ_ant_roots; α1r, α2r, α3r)
       #energy_activation, # E0 = f(T_ant_soil; α4, Estar)
       #source_temperature_coeff, 
       #root_source,
       #microbe_source_moisture_coeff,
       #decomposition_potential,
       #soluble_soil_carbon,
       (x, y, p) -> microbe_source(x, x + FT(10.0), y, y, Csom, Cmic, p),
       co2_diffusivity,
       (x, y, p) -> root_source(x, x + FT(10.0), y, y, Cr, p) 
      ]

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


