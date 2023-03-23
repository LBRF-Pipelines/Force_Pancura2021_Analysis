# -*- coding: utf-8 -*-

import os
import csv


def is_numeric(x):
    try:
        x + 2
        return True
    except TypeError:
        return False


def get_metadata(cfs, framevar):
    vals = []
    as_number = is_numeric(cfs['frames'][0].vars[framevar])
    for i in range(len(cfs['frames'])):
        v = cfs['frames'][i].vars[framevar]
        if as_number:
            v = float(v)
        vals.append(v)
    return vals


def write_metadata(cfs, outpath, framevars=[]):

    filename = cfs['header']['file']['filename']
    frame_num = 0

    header = ['filename', 'frame'] + framevars
    rows = []
    for i in range(len(cfs['frames'])):
        row = [filename, i + 1]
        for var in framevars:
            row.append(cfs['frames'][i].vars[var])
        rows.append(row)

    with open(outpath, 'w+', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def write_channel(cfs, outpath, chan):

    header = ['frame', 'time'] + [chan]

    with open(outpath, 'w+', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for i in range(len(cfs['frames'])):
            frame = cfs['frames'][i]
            time = 0
            dt = frame.sample_rates[chan]
            time_sigfigs = len(str(dt).split('.')[-1])
            for val in frame.data[chan]:
                writer.writerow([i + 1, round(time, time_sigfigs), val])
                time += dt
