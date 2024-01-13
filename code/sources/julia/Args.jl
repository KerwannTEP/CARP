using ArgParse

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--nbAvr"
    help = "Number of sampling points for the radius orbit-averaging integral"
    arg_type = Int64
    default = 20
    "--nbu"
    help = "Number of sampling points for the radius orbit-averaging integral for Jr,L gradients"
    arg_type = Int64
    default = 300
    "--q"
    help = "Anisotropy parameter q for the Plummer model."
    arg_type = Float64
    default = 0.0
    "--a"
    help = "Rotation parameter for the Lynden-Bell Demon."
    arg_type = Float64
    default = 0.0
    "--nbw"
    help = "Number of sampling points for the 3D Rosenbluth integrals: w-integral"
    arg_type = Int64
    default = 20
    "--nbphi"
    help = "Number of sampling points for the 3D Rosenbluth integrals: phi-integral"
    arg_type = Int64
    default = 20
    "--nbvartheta"
    help = "Number of sampling points for the 3D Rosenbluth integrals: vartheta-integral"
    arg_type = Int64
    default = 20

end
parsed_args = parse_args(tabargs)
