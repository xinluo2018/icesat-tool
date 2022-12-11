# author: xin luo
# create: 2021.9.15; modify: 2022.8.12
# des: read in and write out reorganized atl03 data


'''
# des: 1) read atl03 data and select variables to write out; 
#      2) add orbit orientation (ascending/descending) variable.
example:
    python read_atl03.py ./input/path/*.h5 -o /output/path/dir -n 4
to do: 
    filter the atl03 by using quality variables? 
'''

import os
import h5py
import numpy as np
import argparse
from astropy.time import Time


def get_args():
    description = "read ICESat-2 ATL03 data files by groud track and orbit."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "ifiles", metavar="ifiles", type=str, nargs="+",
        help="input files to read (.h5).")
    parser.add_argument(
        '-o', metavar=('outdir'), dest='outdir', type=str, nargs=1,
        help='path to output folder', 
        default=[""])
    parser.add_argument(
        "-n", metavar=("njobs"), dest="njobs", type=int, nargs=1,
        help="number of cores to use for parallel processing", 
        default=[1])
    return parser.parse_args()

def gps2dyr(time):
    """ Converte from GPS time to decimal years. """
    time = Time(time, format="gps")
    time = Time(time, format="decimalyear").value
    return time

def orbit_type(time, lat):
    """
    des: determines ascending and descending tracks
         through testing whether lat increases when time increases.
    arg:
        time, lat: time and latitute of the pohton points.
    return:
        i_asc, i_des: track of the photon points, 1-d data consist of True/Talse. 
    """
    tracks = np.zeros(lat.shape)
    # set track values, !!argmax: the potential turn point of the track
    tracks[0: np.argmax(np.abs(lat))] = 1
    i_asc = np.zeros(tracks.shape, dtype=bool)

    # loop through unique tracks: [0]/[1]/[0,1]
    for track in np.unique(tracks):
        (i_track,) = np.where(track == tracks)
        if len(i_track) < 2:  # number of photon points of specific track i 
            continue
        i_time_min, i_time_max  = time[i_track].argmin(), time[i_track].argmax()
        lat_diff = lat[i_track][i_time_max] - lat[i_track][i_time_min]
        # Determine track type
        if lat_diff > 0:
            i_asc[i_track] = True
    return i_asc, np.invert(i_asc)


def read_atl03(file_in, dir_out):
    ''' 
    des: 
        selecting atl03 variables.
        spot 1,2,3,4,5,6 are always corresponding to beam 1,2,3,4,5,6 
        and spot 1,3,5 are strong beams, spot 2, 4, 6 are weak beams.
        users can add the interested variables by themself.
    arg: 
        file_in: atl03 file, .h5 format
        dir_out: path to save the splitted atl03 data
    return:
        selected atl03 variables
    '''

    # Create dictionary for saving output variables
    out_keys = ['lon', 'lat', 'h', 't_dyr', 'cycle', 'rgt', \
                        'beam_type', 'spot', 'orbit_type']   ## output variables
    d, d_update = {}, {}      
    for key in out_keys: 
        d[key]=np.array([]); 
        d_update[key]=np.array([]);

    #-----------------------------------#
    # 1) Read data for beams         #
    #-----------------------------------#
    group = ["gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"]
    ## loop for groups
    for k in range(len(group)):
        with h5py.File(file_in, "r") as fi:
            try:
                ## group varibales:
                d['lat'] = fi[group[k] + "/heights/lat_ph"][:]
                d['lon'] = fi[group[k] + "/heights/lon_ph"][:]
                d['h'] = fi[group[k] + "/heights/h_ph"][:]
                d['t_dt'] = fi[group[k] + "/heights/delta_time"][:]
                ## dset varibales
                d['tref'] = fi["/ancillary_data/atlas_sdp_gps_epoch"][:]
                d['cycle'] = fi["/orbit_info/cycle_number"][:] * np.ones(len(d['lat'])).astype(np.int8)
                d['rgt'] = fi["/orbit_info/rgt"][:] * np.ones(len(d['lat'])).astype(np.int32)
                d['signal_conf_ph'] = fi[group[k] + "/heights/signal_conf_ph"][:]
                ## group attributes
                beam_type = fi[group[k]].attrs["atlas_beam_type"].decode()
                spot_number = fi[group[k]].attrs["atlas_spot_number"].decode()   # 
            except:
                print(("missing group:", group[k]))
                print(("in file:", file_in))
                continue

        ## set beam type: 1 -> strong, 0 -> weak
        if beam_type == "strong":
            d['beam_type'] = np.ones(d['lat'].shape).astype(np.int8)
        else:
            d['beam_type'] = np.zeros(d['lat'].shape).astype(np.int8)

        #----------------------------------------------------#
        # 3) obtain orbit orientation with time: 
        #    ascending -> 1, descending -> 0                 #
        #----------------------------------------------------#
        ### --- creating array of spot numbers
        d['spot'] = float(spot_number) * np.ones(d['lat'].shape).astype(np.int8)
        t_gps = d['t_dt'] + d['tref']
        d['t_dyr'] = gps2dyr(t_gps)       #  time in decimal years
        ### --- obtain orbit type
        (i_asc, i_des) = orbit_type(d['t_dyr'], d['lat'])    #  track type (asc/des)        
        d['orbit_type'] = np.empty_like(d['lat'], dtype=int)
        d['orbit_type'][i_asc] = 1     # ascending
        d['orbit_type'][i_des] = 0     # descending

        #----------------------------------------------------#
        # 4) Update the variables   #
        #----------------------------------------------------#
        ### ---- Update the variables in dictionary
        for key in out_keys:
            d_update[key] = np.append(d_update[key], d[key])

    #------------------------------------------#
    # 4) Writting out the selected data        #
    #------------------------------------------#
    name, ext = os.path.splitext(os.path.basename(file_in))
    file_out = os.path.join(dir_out, name + "_" + "readout" + ext)
    with h5py.File(file_out, "w") as f_out:
        [f_out.create_dataset(key, data=d_update[key]) for key in out_keys]
    print('written file:', file_out)
    return

if __name__ == '__main__':

    ### ---- read input from command line
    args = get_args()
    ifiles = args.ifiles
    dir_out = args.outdir[0]
    njobs = args.njobs[0]

    if njobs == 1:
        print("running in serial ...")
        [read_atl03(f, dir_out) for f in ifiles]
    else:
        print(("running in parallel (%d jobs) ..." % njobs))
        from joblib import Parallel, delayed
        Parallel(n_jobs=njobs, verbose=5)(delayed(read_atl03)(f, dir_out) for f in ifiles)


