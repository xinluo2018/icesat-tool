## author: Fernando Paolo
## modify: xin luo, 2022.10.7
## des: find and compute crossover values for the icesat2 data. 
## note: !!Due to icesat2 have 6 beams, this script is not suitable for another altimetry data.

import numpy as np
import pyproj
import pandas as pd
import warnings
from scipy import stats
# Ignore all warnings
warnings.filterwarnings("ignore")

def intersect(x_up, y_up, x_down, y_down):
    """
    !!! more fast, but difficult to understand
    reference: 
        https://stackoverflow.com/questions/17928452/
        find-all-intersections-of-xy-data-point-graph-with-numpy
    des: Find orbit crossover locations through solving the equation: 
        p0 + s*(p1-p0) = q0 + t*(q1-q0); p and q are descending and ascending points respectively.
        ---> s*(p1-p0)-t*(q1-q0) = q0-p0
        if s and t belong to [0,1], p and q actually do intersect.
        !! in order to speed up calculation, this code vectorizing solution of the 2x2 linear systems
    input:
        x_down, y_down: coord_x and coord_y of the descending points
        x_up, y_up: coord_x, coord_y of the ascending points.
    retrun:
          coord_x, coord_y: the intersection points. 
          idx_pre_down, idx_pre_up: the index of the pre-point of the down and up tracks, respectively.
    """
    p = np.column_stack((x_down, y_down))   # coords of the descending points
    q = np.column_stack((x_up, y_up))       # coords of the ascending points
    (p0, p1, q0, q1) = p[:-1], p[1:], q[:-1], q[1:]   # remove first/last row respectively
    # (num_uppoints, 2) - (num_dpoints, 1, 2), array broadcast, dim: (num_dpoints, num_uppoints, 2)
    rhs = q0 - p0[:, np.newaxis, :]    
    mat = np.empty((len(p0), len(q0), 2, 2))  # dim: (p_num, q_num, dif((x, y)), orbit(down,up))
    mat[..., 0] = (p1 - p0)[:, np.newaxis]  # dif (x_down,y_down) between point_down and previous point_down
    mat[..., 1] = q0 - q1      #  dif (x_up, y_up) between point_up and previous point_up
    mat_inv = -mat.copy()
    mat_inv[..., 0, 0] = mat[..., 1, 1]   # exchange between x_dif and y_dif, down and up
    mat_inv[..., 1, 1] = mat[..., 0, 0]
    det = mat[..., 0, 0] * mat[..., 1, 1] - mat[..., 0, 1] * mat[..., 1, 0]
    mat_inv /= det[..., np.newaxis, np.newaxis]    # ???
    params = mat_inv @ rhs[..., np.newaxis]        # 
    intersection = np.all((params >= 0) & (params <= 1), axis=(-1, -2))  #
    p0_s = params[intersection, 0, :] * mat[intersection, :, 0]
    xover_coords = p0_s + p0[np.where(intersection)[0]]
    idx_pre_down, idx_pre_up = np.where(intersection)[0], np.where(intersection)[1]
    xover_x, xover_y = xover_coords[:,0], xover_coords[:,1]
    ## remove unreasonable xover points. avoid the crossover point located between 
    ## the one track end point and another track start point
    d_up_mean=np.mean(abs(x_up[1:]-x_up[0:-1]))
    d_down_mean=np.mean(abs(x_down[1:]-x_down[0:-1]))
    x_up_pre, x_up_post = x_up[idx_pre_up], x_up[idx_pre_up+1]
    x_down_pre, x_down_post = x_down[idx_pre_down], x_down[idx_pre_down+1]
    d_up_post_pre, d_down_post_pre = abs(x_up_post-x_up_pre), abs(x_down_post-x_down_pre)
    idx_save = np.argwhere((d_up_post_pre<d_up_mean*3) & (d_down_post_pre<d_down_mean*3)).flatten()
    xover_x, xover_y, idx_pre_up, idx_pre_down = xover_x[idx_save], \
                                    xover_y[idx_save], idx_pre_up[idx_save], idx_pre_down[idx_save]
    return xover_x, xover_y, idx_pre_up, idx_pre_down



def interp_xover(x0, y0, z0, x1, y1, z1, xi, yi):
    '''
    des: interpolation of oxver point by using the pre- and post-points.  
         this function can be used to interpolate the height and time of the oxver point. 
    args:
        x0,y0: the coordinates of pre-point.
        x1,y1: the coordiantes of post-point.
        xi,yi: the coordinates of interpolation point.
        z0: the value of pre-point.
    return:
        zi: the interpolation value.
    ''' 
    ### d: the distance bettwen the pre- and post-points.
    d = np.sqrt(np.square(x1-x0) + np.square(y1-y0))
    di = np.sqrt(np.square(xi-x0) + np.square(yi-y0))
    zi = z0 + (di/d)*(z1-z0)
    return zi


def get_bboxs_id(x, y, xmin, xmax, ymin, ymax, dxy, buff):
    """
    des: get binn box (bbox) id of each points (x,y).
    arg:
        x,y: coordinates of the photon points
        xmin/xmax/ymin/ymax: must be in grid projection: stereographic (m).
        dxy: grid-cell size.
        buff: buffer region, unit is same to x, y
    return:
        the index of each points corresponding to the generated bins.  
    """
    dx, dy = dxy
    # Number of tile edges on each dimension
    Nns = int(np.abs(ymax - ymin) / dy) + 1
    New = int(np.abs(xmax - xmin) / dx) + 1
    # Coord of tile edges for each dimension
    xg = np.linspace(xmin-buff, xmax+buff, New)
    yg = np.linspace(ymin-buff, ymax+buff, Nns)
    # Indicies for each points
    bboxs_id = stats.binned_statistic_2d(x, y, np.ones(x.shape), 'count', bins=[xg, yg]).binnumber
    return bboxs_id

def coor2coor(srs_from, srs_to, x, y):
    """ see utils/transform_xy
    """
    srs_from = pyproj.Proj(int(srs_from))
    srs_to = pyproj.Proj(int(srs_to))
    return pyproj.transform(srs_from, srs_to, x, y, always_xy=True)

def xover_icesat2(lon_as, lat_as, t_as, h_as, spot_as, 
            lon_des, lat_des, t_des, h_des, spot_des, 
            proj, tile_dxy=[20, 20], buff=2):
    """ 
    des: find and compute crossover values. 
    arg:
        lon_as, lat_as: coordinates of ascending points.
        t_as: time of ascending points.
        h_as: height of ascending points
        spot_as: groud track (0-5) of the ascending points.
        lon_des, lat_des: coordinates of descending points.
        t_des: time of descending points.
        h_des: height of descending points
        spot_des: groud track (0-5) of the descending points.
        proj: projection (espg number).
        tile_dxy: width/height of the generated tile. For speeding up processing. unit:km
        buff: buffer of the tile. unit: km
    return:
        out: 
    """
    if proj == "4326":
        raise ValueError("proj can't be 4326")
    tile_dxy[0] *= 1e3
    tile_dxy[1] *= 1e3
    buff *= 1e3
    ######## -------- 1. find the xover points -------- #####
    # Transform to wanted coordinate system
    (x_as, y_as) = coor2coor(4326, proj, lon_as, lat_as)
    (x_des, y_des) = coor2coor(4326, proj, lon_des, lat_des)

    # spatial range (m)
    xmin = max(np.nanmin(x_as), np.nanmin(x_des))
    xmax = min(np.nanmax(x_as), np.nanmax(x_des))
    ymin = max(np.nanmin(y_as), np.nanmin(y_des))
    ymax = min(np.nanmax(y_as), np.nanmax(y_des))

    print('tiling asc/des data...') 
    # get binned boxes index (bin indices corresponding to each (x,y))
    # here the bin is tile.
    id_bboxs_as = get_bboxs_id(x_as, y_as, xmin, xmax, ymin, ymax, tile_dxy, buff)     # box id of the ascending photon
    id_bboxs_des = get_bboxs_id(x_des, y_des, xmin, xmax, ymin, ymax, tile_dxy, buff)  # box id of the descending photon
    # Copy box for convenience
    id_bboxs = id_bboxs_as
    # Initiate output container (dictionary)
    out = []    
    ibox = np.unique(id_bboxs)
    print('searching ibox (sub-tile):', ibox)    
    num_box = 0              #  counter: count of the bins valid

    print('computing crossovers ...')
    #######   for bin in bins:
    #######       for track_as in tracks_as:    ## track is beam of the icesat2.
    #######           for track_des in tracks_des: 
    #######               find the xover_points.
    # loop through each unique bin (tile)
    # k is the unique bin index.
    for k in ibox:   
        ibox_as, = np.where(id_bboxs_as == k)    # ibox_as is the box index of data points 
        ibox_des, = np.where(id_bboxs_des == k)
        # Extract points in the bin
        # ascending orbit
        spot_as_ibox = spot_as[ibox_as]
        x_as_ibox = x_as[ibox_as]
        y_as_ibox = y_as[ibox_as]
        h_as_ibox = h_as[ibox_as]
        t_as_ibox = t_as[ibox_as]
        # descending orbit        
        spot_des_ibox = spot_des[ibox_des]   
        x_des_ibox = x_des[ibox_des]
        y_des_ibox = y_des[ibox_des]
        h_des_ibox = h_des[ibox_des]
        t_des_ibox = t_des[ibox_des]

        # get unique tracks
        spot_as_ibox_unique = np.unique(spot_as_ibox)
        spot_des_ibox_unique = np.unique(spot_des_ibox)

        # Test if tile has no data 
        # cause bboxs = bboxs1, len(orbits1) !=1, len(orbits2) could be 0.
        if len(spot_as_ibox) == 0 or len(spot_des_ibox) == 0:
            continue

        # ---- loop through tracks (ground track in the specific bin) 
        # --> ascending tracks
        for ispot_as_ibox in spot_as_ibox_unique:
            ## i_trk_: point index of the specific groud track.
            i_as_spot = spot_as_ibox == ispot_as_ibox 
            ## extract points from the specific orbit (a specific track)
            x_as_ispot = x_as_ibox[i_as_spot]
            y_as_ispot = y_as_ibox[i_as_spot]
            t_as_ispot = t_as_ibox[i_as_spot]
            h_as_ispot = h_as_ibox[i_as_spot]
            
            # Loop through groud tracks (1-6) in specific bin 
            # --> descending tracks
            for ispot_des_ibox in spot_des_ibox_unique:
                # index of data points of specific track in file 2
                i_des_spot = spot_des_ibox == ispot_des_ibox
                # extract points from a specific orbit (groud track)
                x_des_ispot = x_des_ibox[i_des_spot]
                y_des_ispot = y_des_ibox[i_des_spot]
                t_des_ispot = t_des_ibox[i_des_spot]
                h_des_ispot = h_des_ibox[i_des_spot]                
                # Test length of vector
                if len(x_as_ispot) < 3 or len(x_des_ispot) < 3: 
                    continue

                # exact crossing points between two tracks of ascending/descending files.
                xi, yi, idx_pre_as, idx_pre_des = intersect(x_as_ispot, y_as_ispot, x_des_ispot, y_des_ispot)
                # ensure the xover points exit
                if len(xi) == 0: continue
                
                ## time of the interpolation photon.
                ti_as = interp_xover(x0=x_as_ispot[idx_pre_as], y0=y_as_ispot[idx_pre_as], z0=t_as_ispot[idx_pre_as], 
                                     x1=x_as_ispot[idx_pre_as+1], y1=y_as_ispot[idx_pre_as+1], z1=t_as_ispot[idx_pre_as+1], 
                                     xi=xi, yi=yi)
                ti_des = interp_xover(x0=x_des_ispot[idx_pre_des], y0=y_des_ispot[idx_pre_des], z0=t_des_ispot[idx_pre_des], 
                                      x1=x_des_ispot[idx_pre_des+1], y1=y_des_ispot[idx_pre_des+1], z1=t_des_ispot[idx_pre_des+1], 
                                      xi=xi, yi=yi)

                ## height of the interpolation photon.
                hi_as = interp_xover(x0=x_as_ispot[idx_pre_as], y0=y_as_ispot[idx_pre_as], \
                                z0=h_as_ispot[idx_pre_as], x1=x_as_ispot[idx_pre_as+1], \
                                y1=y_as_ispot[idx_pre_as+1], z1=h_as_ispot[idx_pre_as+1], xi=xi, yi=yi)
                hi_des = interp_xover(x0=x_des_ispot[idx_pre_des], y0=y_des_ispot[idx_pre_des], \
                                z0=h_des_ispot[idx_pre_des], x1=x_des_ispot[idx_pre_des+1], \
                                y1=y_des_ispot[idx_pre_des+1], z1=h_des_ispot[idx_pre_des+1], xi=xi, yi=yi)

                # Create output array
                out_i = np.full((10, len(xi)), np.nan)
                # Compute differences and save parameters
                out_i[0]  = xi              # crossover points coord_x
                out_i[1]  = yi              # ... coord_y
                out_i[2]  = hi_as           # interpolated height by ascending track
                out_i[3]  = hi_des          # interpolated height by descending track
                out_i[4]  = ti_as           # interpolated time by ascending track
                out_i[5]  = ti_des           # interpolated time by descending track
                out_i[6] = ispot_as_ibox      # groud track of ascending file
                out_i[7] = ispot_des_ibox     # groud track of decending file
                out_i[8]  = hi_as - hi_des    ## height difference between ascending and descending interpolations
                out_i[9]  = ti_as - ti_des    ## time difference between ...        
                # Add to list
                out.append(out_i)
        num_box += 1
    # Change back to numpy array
    # Test if output container is empty 
    if len(out) == 0:
        print('no crossovers found!')
        return 
    out = np.concatenate(out, axis=1).T

    # Remove invalid rows
    out = out[~np.isnan(out[:,2]), :]     # out[:,2]: height difference
    # Test if output container is empty 
    if len(out) == 0:
        print('no crossovers found!')
        return 
    print('number of crossovers:', len(out))
    # remove the two id columns if they are empty, out[:,-1]: orb_id2，out[:,-2]: orb_id1
    out = out[:,:-2] if np.isnan(out[:,-1]).all() else out
    # Transform coords back to lat/lon
    out[:,0], out[:,1] = coor2coor(proj, '4326', out[:,0], out[:,1])
    out_df = pd.DataFrame(out, columns=['o_lon', 'o_lat', 'oh_as', 'oh_des', 'ot_as', 'ot_des', \
                                    'ospot_as', 'ospot_des', 'oh_dif', 'ot_dif'])
    return out_df

