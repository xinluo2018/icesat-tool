"""
author: Fernando Paolo, 
modify: xin luo, 2022.8.2.  
des: 
  read and write out icesat1-gla14 data.
  main processing: 1) height corrections, 2) selected valid points, 
                   3) selected variables. 3) add oribt type
example:
  python read_glah14.py ./input/path/*.h5 -o /output/path/dir -n 4 
"""

import os
import h5py
import numpy as np
import argparse

def utc2dyr(time):
  # compute correct time: add back 'year 2000 + 12 hours' in secs
  # convert time to decimal years
  time += (2000 * 365.25 * 24 * 3600.) + (12 * 3600.)
  time = time/(365.25*24*3600.)
  return time

def orbit_type(time, lat):
    """
        Determines ascending and descending tracks.
        Defines unique tracks as segments with time breaks > tmax,
        and tests whether lat increases or decreases w/time.
    """    
    tracks = np.zeros(lat.shape)     # Generate track segment
    tracks[0:np.argmax(np.abs(lat))] = 1   # Set values for segment
    i_asc = np.zeros(tracks.shape, dtype=bool)  # Output index array
    # Loop trough individual tracks
    for track in np.unique(tracks):
        i_track, = np.where(track == tracks)  # Get all points from an individual track
        # Test tracks length
        if len(i_track) < 2:
            continue
        # Test if lat increases (asc) or decreases (des) w/time
        i_min = time[i_track].argmin()
        i_max = time[i_track].argmax()
        lat_diff = lat[i_track][i_max] - lat[i_track][i_min]        
        # Determine track type
        if lat_diff > 0:
            i_asc[i_track] = True
    # Output index vector's
    return i_asc, np.invert(i_asc)
    
def get_args():
    description = "read ICESat-2 ATL06 data files by groud track and orbit."
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

def read_glah14(file_in, dir_out):
    ''' 
    des: 
        elevation corrections, add oribt type;
    arg: 
        file_in: glas14 file, .H5 format
        path_out: path to save the splitted glas14 data
    return:
        the corrected glas14 height and the selected glas14 data. 
    '''


    d = {}  # Dictionary for saving selected variables
    #------------------------------------------#
    # 1) read in the 40HZ data        #
    #------------------------------------------#
    with h5py.File(file_in, 'r') as f_gla14:

        d['lat'] = f_gla14['Data_40HZ/Geolocation/d_lat'][:]
        d['lon'] = f_gla14['Data_40HZ/Geolocation/d_lon'][:]
        d['t_utc'] = f_gla14['Data_40HZ/Time/d_UTCTime_40'][:]  # secs since 2000-01-01 12:00:00 UTC
        d['h'] = f_gla14['Data_40HZ/Elevation_Surfaces/d_elev'][:] 
        d['h_flg'] = f_gla14['Data_40HZ/Quality/elev_use_flg'][:]
        d['sat_flg'] = f_gla14['Data_40HZ/Quality/sat_corr_flg'][:]
        d['att_flg'] = f_gla14['Data_40HZ/Quality/sigma_att_flg'][:]
        d['h_sat'] = f_gla14['Data_40HZ/Elevation_Corrections/d_satElevCorr'][:]  # saturation cor [m]
        d['num_pk'] = f_gla14['Data_40HZ/Waveform/i_numPk'][:]
        d['h_ellip'] = f_gla14['Data_40HZ/Geophysical/d_deltaEllip'][:]
        d['h_tide'] = f_gla14['Data_40HZ/Geophysical/d_ocElv'][:]  # ocean tide [m]
        d['h_load'] = f_gla14['Data_40HZ/Geophysical/d_ldElv'][:]  # load tide [m]  
        track_01Hz = f_gla14['Data_1HZ/Geolocation/i_track'][:]   # 1Hz Track 
        d['track'] = track_01Hz.repeat(40)  # Get 40 Hz tracks

        d['t_dyr'] = utc2dyr(d['t_utc'])

        #------------------------------------------------------------#
        # 2) selected the valid data and perform height corrections #
        #------------------------------------------------------------#

        idx, = np.where(
            (np.abs(d['h']) < 1e10) &
            (np.abs(d['h_sat']) < 1e10) &
            (np.abs(d['lat']) <= 90) &
            (np.abs(d['lon']) <= 3600) &
            (d['h_flg'] == 0) &
            (d['sat_flg'] <= 2) &
            (d['att_flg'] == 0) &
            (d['num_pk'] == 1))

        # Check if no valid pts
        if len(idx) == 0:
            print(('no valid pts:', file_in))
            return

        ## Selected the valid points.
        for k in list(d.keys()):   
            d[k] = d[k][idx]

        ## Corrections
        d['h'] += d['h_tide'] + d['h_load']  # apply tides cor
        d['h'] += d['h_sat']      # Apply saturation cor
        d['h'] -= d['h_ellip']    # Convert ellipsoid: h_TP -> h_WGS84


        #------------------------------------------#
        # 3) add orbit type         #
        #------------------------------------------#
        name, ext = os.path.splitext(os.path.basename(file_in))
        (i_asc, i_des) = orbit_type(d['t_dyr'], d['lat'])
        d['orbit'] = np.empty_like(d['lat'], dtype=int)
        d['orbit'][i_asc] = 1  # 1 -> ascending
        d['orbit'][i_des] = 0  # 0 -> descending

        #------------------------------------------#
        # 4) Writting out the selected data        #
        #------------------------------------------#
        out_keys = ['lon', 'lat', 'h', 't_dyr', 'track', 'orbit']   ## output variables
        file_out = os.path.join(dir_out, name + '_readout' + '.h5')   ## note: change the extention name .H5 to .h5.
        with h5py.File(file_out, "w") as f_out:
            [f_out.create_dataset(key, data=d[key]) for key in out_keys]
        print('written file:', (file_out))

if __name__ == '__main__':
    ### ---- read input from command line
    args = get_args()
    ifiles = args.ifiles
    dir_out = args.outdir[0]
    njobs = args.njobs[0]

    if njobs == 1:
        print("running in serial ...")
        [read_glah14(f, dir_out) for f in ifiles]
    else:
        print(("running in parallel (%d jobs) ..." % njobs))
        from joblib import Parallel, delayed
        Parallel(n_jobs=njobs, verbose=5)(delayed(read_glah14)(f, dir_out) for f in ifiles)

