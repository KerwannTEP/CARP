using ArgParse

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--parallel"
    help = "Parallel computation: yes/no"
    arg_type = String
    default = "yes"
    "--nbK"
    help = "Number of sampling points for the 3D Rosenbluth integrals"
    arg_type = Int64
    default = 100
    "--nbAvr"
    help = "Number of sampling points for the radius orbit-averaging integral"
    arg_type = Int64
    default = 100
    "--nbAvrTh"
    help = "Number of sampling points for the theta orbit-averaging integral"
    arg_type = Int64
    default = 100
    "--nbu"
    help = "Number of sampling points for the radius orbit-averaging integral for Jr,L gradients"
    arg_type = Int64
    default = 300
    "--q"
    help = "Anisotropy parameter q for the Plummer model."
    arg_type = Float64
    default = 1.0
    "--a"
    help = "Rotation parameter for the Lynden-Bell Demon."
    arg_type = Float64
    default = 0.1
    "--g"
    help = "Rotation fonction g(cos I)"
    arg_type = String
    default = "SIGN"
    "--a_err"
    help = "a parameter for err-g function"
    arg_type = Float64
    default = 1.0
    "--nbw"
    help = "Number of sampling points for the 3D Rosenbluth integrals: w-integral"
    arg_type = Int64
    default = 100
    "--nbphi"
    help = "Number of sampling points for the 3D Rosenbluth integrals: phi-integral"
    arg_type = Int64
    default = 100
    "--nbvarphi"
    help = "Number of sampling points for the 3D Rosenbluth integrals: varphi-integral"
    arg_type = Int64
    default = 100

end
parsed_args = parse_args(tabargs)
