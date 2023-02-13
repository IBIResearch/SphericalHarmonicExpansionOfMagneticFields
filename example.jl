##################################
# Selection Field in MPI (2 T/m) #
##################################
# Example for the paper
# "Unique Compact Representation 
#  of Magnetic Fields  
#  using Truncated Solid 
#  Harmonic Expansions"
# by M. Boberg, T. Knopp, and
#    M. Möddel
##################################
# Julia: 1.8

#########
# Setup #
#########
# Activate local environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Packages
using SphericalHarmonicExpansions, HDF5
using NLsolve # used to find FFPs
using PyPlot # used for visualization
using LazyArtifacts # used to download the data

# Variabels for polynomial expansions
@polyvar x y z

# Close all previous figures
close("all")

# include utils
include("utils/magneticField.jl") # magneticField() to calculate coefficients of magnetic field
include("utils/findFFP.jl") # findFFP() to find FFPs of magnetic field
include("utils/errorPropagation.jl") # errorPropagation() to propagate the measurement error to the coefficients
include("utils/plotMagneticField.jl") # plotMagneticField() for visualization of the magnetic field
include("utils/plotCoeffs.jl") # plotCoeffs() to visualize the coefficients

################
# 1. Load Data #
################
# path to file
pathFile = "data/Gradient_2Tpm.h5"

# read data from hdf5 file
field, fieldsError, coords, R, center = h5open(pathFile, "r") do file
    field = read(file,"/fields") 		# measured field
    fieldsError = read(file,"/fieldsError") 	# field error stemming from the gaussmeter
    coords = read(file,"positions") 		# measured positions (shifted and scaled t-design)
    R = read(file,"positionsTDesignRadius")	# radius of the measured ball
    center = read(file,"positionsCenter")	# center of the measured ball
    return field, fieldsError, coords, R, center
end


###########################
# 2. Initial Coefficients #
###########################
coeffs, expansion, func = magneticField(coords,field, R,center,4,x,y,z)


######################
# 3. Post Processing #
######################
## 3.1 Translation in the Hall sensor
# Translation vectors for all three directions
v = [-0.0018   -0.0018   -0.0018;
      0.0       0.00208   0.0;
     -0.00208   0.0       0.0]
# Translation
for c=1:size(coeffs,2), j = 1:3
  coeffs[j,c] = SphericalHarmonicExpansions.translation(coeffs[j,c],v[:,j])
  expansion[j,c] = sphericalHarmonicsExpansion(coeffs[j,c],x,y,z);
end

## 3.2  Translation to FFP 
# Finding the FFP
ffp = findFFP(expansion, x, y, z)
# Translation (each coefficient by ffp[:,c])
for c=1:size(coeffs,2), j = 1:3
  coeffs[j,c] = SphericalHarmonicExpansions.translation(coeffs[j,c],ffp[:,c])
  expansion[j,c] = sphericalHarmonicsExpansion(coeffs[j,c],x,y,z);
end

## 3.3 Error Propagation
coeffs, coeffsError, expansion, func = errorPropagation(coeffs,coords,R,center,fieldsError, x, y, z, ffp=ffp)

###############
# 4. Plotting #
###############
# 4.1 Plotting of the magnetic field
# Plotting three orthogonal planes through the FFP with marked measured sphere
fignumber = plotMagneticField(func[:,1], R, center=ffp[:,1])

# 4.2 Plotting of the coefficients
# Plot of the unnormalized coefficients up to degree 4
fignumber = plotCoeffs(coeffs, 4, fignumber=fignumber, separate=false, fontsize=14, R=R)
# Plot of the normalized coefficientc up to degree 2 with propagated error
fignumber = plotCoeffs(coeffs, 2, fignumber=fignumber, separate=true, plotError=true, coeffsError=coeffsError)
