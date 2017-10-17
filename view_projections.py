import genfire.fileio
import matplotlib.pyplot as plt
import sys

# set default projection file if none was provided
if len(sys.argv) < 2:
	file_projections = 'data/projections_FePt.mat'
else:
	# verify the file exists
	from os.path import isfile
	file_projections = sys.argv[1]
	if not isfile(file_projections):
		raise IOError("File {} does not exist!".format(file_projections))

arr = genfire.fileio.readVolume(file_projections)

for i in range(arr.shape[2]):
    plt.figure(42)
    plt.imshow(arr[:, :, i])
    plt.title("Projection #{}/{}\nPress any key to advance or click the mouse to exit".format(i, arr.shape[2]))
    plt.draw()
    plt.pause(1e-6) # forces figure to be rendered
    if not plt.waitforbuttonpress(None):
    	exit(1)
