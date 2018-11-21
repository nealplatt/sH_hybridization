
cd C:\Users\nplatt\Desktop\2beb168a6f53090f

reinitialize
set matrix_mode, 1
set movie_panel, 1

# General settings
hide all
set valence, 1
valence guess, all
set specular, off
set depth_cue, 1
bg_color white
space rgb


# Load data
load final.casp.pdb

# General settings
hide all
set valence, 1
valence guess, all
set specular, off
set depth_cue, 1
bg_color white
space rgb


show ribbon

color gray30


#highlight active sites
sites = ["resi 240", "resi 241", "resi 242", "resi 243", "resi 244"] 
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "blue", i + " and name ca")


#highlight pocket sites
sites = ["resi 396", "resi 397", "resi 398", "resi 399", "resi 400", "resi 401", "resi 405", "resi 415", "resi 423", "resi 424", \
"resi 425", "resi 426", "resi 427", "resi 428", "resi 429", "resi 430", "resi 433", "resi 442", "resi 443", "resi 444", \
"resi 445", "resi 446", "resi 447", "resi 448", "resi 450", "resi 455", "resi 456", "resi 515", "resi 557", "resi 558", \
"resi 559", "resi 560", "resi 561", "resi 563", "resi 564", "resi 565", "resi 566", "resi 567", "resi 568", "resi 569", \ 
"resi 576", "resi 577", "resi 578", "resi 579", "resi 612", "resi 631", "resi 632", "resi 633", "resi 634", "resi 637"]
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_scale", "0.5") \
	cmd.set("sphere_color", "green", i + " and name ca")

#highlight mutated sites
sites = ["resi 94", "resi 95", "resi 107", "resi 169", "resi 296", "resi 309", "resi 565", "resi 595", "resi 710", "resi 711"] 
for i in sites: \
	cmd.show("sphere", i + " and name ca") \
	cmd.set("sphere_color", "red", i + " and name ca")
