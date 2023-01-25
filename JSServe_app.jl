using JSServe
using parameterization_viz 

my_app = param_dashboard(
                         SoilCO2ModelParameters,
                         Dict("CO2 production" => (d1, d2, p) -> microbe_source(d1, d2, 5.0, p), "CO2 diffusivity" => co2_diffusivity),
                         ["T_soil", "M_soil"],
                         ([273, 303], [0.0, 0.5])
                        )

app_name = get(ENV, "HEROKU_APP_NAME", "param-viz")
url = "https://$(app_name).herokuapp.com"
server = JSServe.Server(my_app, "0.0.0.0", parse(Int, ENV["PORT"]))
server.proxy_url = url

wait(server)
