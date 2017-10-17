# Basic script for converting projections from one file type to another
import genfire.fileio

filename_in = 'data/projections_FePt.mat'
filename_out = 'data/projections_FePt.npy'

filename_in = 'data/projections_FePt_noisy.mat'
filename_out = 'data/projections_FePt_noisy.npy'

filename_in = 'data/projections_FePt_misaligned.mat'
filename_out = 'data/projections_FePt_misaligned.npy'

filename_in = 'data/projections_FePt_mw.mat'
filename_out = 'data/projections_FePt_mw.npy'

arr = genfire.fileio.readVolume(filename_in)
genfire.fileio.writeVolume(filename_out, arr)
