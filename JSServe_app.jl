using Pkg
using Dates
Pkg.gc(;collect_delay=Dates.Day(0))

using JSServe
using parameterization_viz 

my_app = param_dashboard(
                         SoilCO2ModelParameters,
                         Dict("CO2 production" => (d1, d2, p) -> microbe_source(d1, d2, 5.0, p), "CO2 diffusivity" => co2_diffusivity),
                         ["T_soil", "M_soil"],
                         ([273, 303], [0.0, 0.5])
                        )

# my_app = App(DOM.div("hello world"))
server = JSServe.Server(my_app, "0.0.0.0", parse(Int, ENV["PORT"]))
#route!(server, "/" => App(DOM.div("hello world")))

wait(server)
