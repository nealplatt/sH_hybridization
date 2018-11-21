# Created in: 2015-11-20


#===========#
# Load data #
#===========#

# setup PyMOL for the movie
reinitialize
set matrix_mode, 1
set movie_panel, 1

# !!! Comment following for preview !!!
set ray_trace_frames=1
set cache_frame=0
mclear

# Load data
load 4MUB.pdb



#===================#
# Redesign molecule #
#===================#

# General settings
hide all
set valence, 1
valence guess, all
set specular, off
set depth_cue, 1
bg_color white
space rgb

# Redesign structure
show ribbon
show sticks, resi 301
show sticks, resi 302

# Redesign color
color gray30, resi -1-109
#color gray70, resi 110-257
color gray30, resi 110-257
color dirtyviolet, resi 301
color magenta, resi 302


# Deleterious mutations
sites = ["resi 35", "resi 142", "resi 119", "resi 225"] 	# The question remains for 225
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "red", i + " and name ca")

# Suspected deleterious mutations
sites = ["resi 171", "resi 179", "resi 183", "resi 225"]
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "orange", i + " and name ca")

# No effect mutations
sites = ["resi 106", "resi 160", "resi 206", "resi 256"]
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "blue", i + " and name ca")

# Suspected no effect mutations
sites = ["resi 3", "resi 25", "resi 44", "resi 50", "resi 171", "resi 57", "resi 66", "resi 71", "resi 74", "resi 152", "resi 174", "resi 176", "resi 180", "resi 188", "resi 189", "resi 204", "resi 203", "resi 223", "resi 227", "resi 243", "resi 244"]
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "marine", i + " and name ca")



#================#
# Rotation movie #
#================#
#run mpng_advanced.py
#mpngAdvanced(800, 800, "ssadh")
#-------------#
# Setup views #
#-------------#

### cut below here and paste into script ###
set_view (\
    -0.451888919,   -0.541844845,    0.708662927,\
     0.857649803,   -0.045342222,    0.512227893,\
    -0.245418206,    0.839253068,    0.485202074,\
     0.000352006,   -0.000050582, -162.170410156,\
   129.527542114,   -0.831031799,   17.474643707,\
    74.684745789,  249.925598145,  -20.000000000 )
### cut above here and paste into script ###
#### cut below here and paste into script ###
#set_view (\
    #-0.451888919,   -0.541844845,    0.708662927,\
     #0.857649803,   -0.045342222,    0.512227893,\
    #-0.245418206,    0.839253068,    0.485202074,\
     #0.000000000,    0.000000000, -162.170410156,\
   #118.462799072,    6.357128143,   12.608121872,\
    #74.549987793,  249.790832520,  -20.000000000 )
#### cut above here and paste into script ###
center

#----------------------#
# Make the movie views #
#----------------------#
mset 1 x180


##---- Python code ----##
python

import os
import imageio	# For processing image into gif
import shutil

if not os.path.exists("test"):
    os.makedirs("test")

myframe = 90
myangle = 360/myframe

for x in range(myframe):
    cmd.frame(x+1)
    cmd.turn("y", myangle)
    cmd.ray(350, 280)		# SD
	#cmd.ray(1024, 576)		# SD
    #cmd.ray(1920,1080)		# HD
    #cmd.png("test/mov_around" + str(x+1) + ".png")
    cmd.png("%s%03d.png" % ("test/mov_around_", x + 1))


filenames = sorted((fn for fn in os.listdir('test/') if fn.endswith('.png')))
images = []
for filename in filenames:
    images.append(imageio.imread("test/" + filename))
imageio.mimsave('movie.gif', images)

shutil.rmtree("test/", ignore_errors=True)

python end
##---- Python code ----##
