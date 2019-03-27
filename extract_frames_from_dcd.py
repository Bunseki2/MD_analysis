#!/usr/bin/env python

# Author: Guillaume Bouvier -- guillaume.bouvier@pasteur.fr
# https://research.pasteur.fr/en/member/guillaume-bouvier/
# 2015-10-27 16:29:46 (UTC+0100)

import MDAnalysis
import numpy
import sys
import os

# This two functions are adapted from:
# http://stackoverflow.com/a/2497565/1679629

def sort_array(initial_array):
    aa = numpy.argsort(initial_array)
    return initial_array[aa]

def unsort_array(sorted_array, initial_array):
    aa = numpy.argsort(initial_array)
    aaa = numpy.argsort(aa)
    if type(sorted_array) is type([]):
        out = []
        for i in aaa:
            out.append(sorted_array[i])
        return out
    else:
        return sorted_array[aaa]

def write_selected_frames(PDB, DCD, frame_list, OUT, selection='all',
                          verbose=False):
    frame_list = numpy.asarray(frame_list)
    if verbose:
        print( "Frame list:\n%s"%frame_list)
    sorted_frame_list = list(sort_array(frame_list)[::-1])
    if verbose:
        print (sorted_frame_list)

    u = MDAnalysis.Universe(PDB, DCD)

    system = u.select_atoms(selection)
    if selection != 'all':
        # Write a pdb of the system that can be used as a topology file to open the
        # resulting trajectory in VMD for example.
        system.write('%s.pdb'%os.path.splitext(OUT)[0])
    timesteps = []
    for ts in u.trajectory:
        try:
            while ts.frame == sorted_frame_list[-1]:
                timesteps.append(system.ts.copy())
                frame_id = sorted_frame_list.pop()
                if verbose:
                    print ("get frame: %d"%frame_id)
        except IndexError:
            break

    timesteps = unsort_array(timesteps, frame_list)
    with MDAnalysis.Writer(OUT, system.n_atoms) as W:
        for ts in timesteps:
            if verbose:
                print ("writing frame: %d"%ts.frame)
            W.write(ts)

if __name__ == '__main__':
    PDB = sys.argv[1]
    DCD = sys.argv[2]
    OUT = sys.argv[3]
    try:
        frame_list_file = sys.argv[4] # ASCII file containing the list of frames.
                                      # One frame id per line and frame ids starts from 0.
                                      # Repetitions of frame are allowed and will be
                                      # repeated in output. The order of the frames will
                                      # be kept in the output.
        frame_list = numpy.genfromtxt(frame_list_file, dtype=int)
    except IndexError: # Read from standard input
        frame_list = numpy.genfromtxt(sys.stdin, dtype=int)
write_selected_frames(PDB, DCD, frame_list, OUT, verbose=True)
