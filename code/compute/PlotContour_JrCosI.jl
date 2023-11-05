using HDF5 # To have access to .hf5 files
using Plots # To be able to plot data
using ArgParse
using LaTeXStrings




##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--q"
    help = "Anisotropy parameter q for the Plummer model."
    arg_type = Float64
    default = 1.0
    "--a"
    help = "Rotation parameter for the Lynden-Bell Demon."
    arg_type = Float64
    default = 0.1

end
parsed_args = parse_args(tabargs)

const qCalc = parsed_args["q"]
const alphaRot = parsed_args["a"]

########################################
# Function that read the .hf5 file
########################################

namefile = "../data/Dump_dFdt_JrCosI_Map_q_"*string(qCalc)*"_a_"*string(alphaRot)*".hf5"

function openData(namefile)
    file = h5open(namefile, "r")
    tabJr = read(file,"tabJr")
    tabCosI = read(file,"tabCosI")
    tabdFdt = read(file,"tabdFdt")

    tabdFdt = tabdFdt .* 10^5

    close(file)
    return tabJr, tabCosI, tabdFdt
end

########################################
# Getting the data
########################################
println("Recovering plot data...")
const tabJr, tabCosI, tabdFdt = openData(namefile)

const maxdata = maximum(abs.(tabdFdt))
const pref = 1.0

########################################
# Plotting the data
########################################
println("Plotting the data...")
p = Plots.contourf(tabCosI,tabJr,tabdFdt,color=:bluesreds, clims= (-pref, pref).* maxdata, title=L" \partial F/\partial t \ [ \ \!\!\!\!\! \times\!\! 10^5]")
Plots.savefig(p,"../graphs/Map_dFdt_JrCosI_q_"*string(qCalc)*"_a_"*string(alphaRot)*".png") # Saves the figure
Plots.display(p)
readline()