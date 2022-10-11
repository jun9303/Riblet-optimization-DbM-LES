import sys, math, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.affinity import scale, translate
from shapely.ops import unary_union
from shapely.validation import make_valid

# local packages
from morph_prelim import lib_gridfunc, lib_gridvalid, lib_griddebug
from morph_prelim import lib_ibmpre, lib_setgrid, lib_ibm_body

Re_tau = 180 # equivalent Re_m   = 5730 when Re_tau = 180
# in the current optimization we assume the above flow conditions 
# Re_tau = u_tau * channelHalfHeight / nu = 180 (given the flat parallel infinite channel)
# this is equivalent to  Re_m = u_bulk * (2*channelHalfHeight) / nu = 5730 (fully turbulent guaranteed)
# ref. Tsukahara et al. (2005) DNS of turbulent channel flow at very low Reynolds numbers
# note: y+ = y_dim.less * Re_tau
# all length inputs are assumed to be in the + scale
N_x = 17  # streamwise direction grid points
N_y = 109  # normal direction grid points
N_z = 321 # spanwise direction grid points
L_x = np.pi # streamwise direction domain length (ref. Choi et al. (1994) Direct numerical simulation of turbulent flow over riblets)
L_y = 2.0 # normal direction domain length (bc scaled by channelHalfHeight). will be adjusted in les_grid_pre
L_z = 1.0 # will be adjusted in les_grid_pre

def arc_len_param_2d(colloc_num, x_coord, y_coord):
    if len(list(x_coord)) != len(list(y_coord)):
        raise TypeError('(arc_len_param_2d) x and y coordinates must have the same length')

    xdiff = np.diff(x_coord)
    ydiff = np.diff(y_coord)
    aleng = np.add.accumulate(np.concatenate((np.array([0.]), np.sqrt(xdiff*xdiff+ydiff*ydiff))))
    alf  = interp1d(aleng, np.r_[1:aleng.size+1], kind='linear')
    div = np.linspace(min(aleng), max(aleng), num=colloc_num, endpoint=True)
    fx = interp1d(np.r_[1:len(list(x_coord))+1], x_coord, kind='linear')
    fy = interp1d(np.r_[1:len(list(y_coord))+1], y_coord, kind='linear')
    x2 = fx(alf(div))
    y2 = fy(alf(div))

    return x2, y2

def coord_save(x_coord, y_coord, filePath, description=''):
    if len(list(x_coord)) != len(list(y_coord)):
        raise TypeError('(coord_save) x and y coordinates must have the same length')

    coord = np.zeros((2, len(list(x_coord))))
    coord[:][0] = x_coord; coord[:][1] = y_coord;
    coord = coord.T

    with open(filePath, 'w') as f:
        f.write('# '+description+'\n')
        np.savetxt(f, coord, delimiter=',', fmt='%15.12f')

    return None

def rib_dbm(baselineDir, weights):
    if baselineDir[-1] != '/':
        baselineDir += '/'
    baselinefiles = glob.glob(baselineDir+'*.txt')
    baselinefiles.sort(key=lambda f: int(re.sub('\D', '', f)))

    if len(list(baselinefiles)) == 0:
        raise RuntimeError('(rib_dbm) no baseline file found in the given filepath')

    if len(list(weights)) > len(list(baselinefiles)):
        raise TypeError('(rib_dbm) too many weight inputs. check the number of baseline')

    baselinecoordinates=[]
    
    for baselinefile in baselinefiles:
        baselinecoordinates.append(np.genfromtxt(baselinefile, delimiter=',', dtype='float64'))

    weights = np.array(weights)

    if np.sum(weights) == 0:
        # raise ValueError('(rib_dbm)sum of the weights must not be zero')
        return [Polygon([(0,0), (1,0), (0,0)]), weights]

    weights_norm = weights / np.sqrt(np.sum(weights*weights))
    if np.sum(weights_norm) < 0:
        weights_norm += 2./len(list(weights))*np.sum(weights_norm)
    # now all weight vector should be on a unit n-sphere over the plane described by sum_w_i = 0
    # the weight surface w_i > 0 does interpolation, while the other surface does extrapolation


    morph = np.zeros(baselinecoordinates[0].shape)

    for wnum, weight in enumerate(weights_norm):
        morph += weight * baselinecoordinates[wnum]

    morph /= np.sum(weights)
    morph.T[1] /= np.max(morph.T[1]) - np.min(morph.T[1])
    # geometric normalization performed C(l(0)) = (0,0), C(l(end)) = (1,0)

    morph.T[1] = np.abs(morph.T[1])
    # negative y treatment -- here we take the strategy of getting abs(y) of the curve

    polygon = Polygon(np.c_[morph.T[0],morph.T[1]])
    # intersection control begins -- create a (multi-)polygon from the morphed curve
    # here we have the outermost boundary of the curve as a resulting shape

    list_parts=[]
    eps = 1e10
    if polygon.geom_type == 'MultiPolygon':
        for pp in polygon.geoms:
            list_interiors=[]
            for interior in pp.interiors:
                p = Polygon(interior)
                
                if p.area > eps:
                    list_interiors.append(interior)
            temp_pol = Polygon(pp.exterior.coords, holes=list_interiors)
            list_parts.append(temp_pol)
        polygon=MultiPolygon(list_parts)
    else:
        list_interiors=[]
        for interior in polygon.interiors:
            p = Polygon(interior)
            if p.area > eps:
                list_interiors.append(interior)
        polygon=Polygon(polygon.exterior.coords, holes=list_interiors)
    # remove any unnecessary inner holes

    polygon = polygon.buffer(1e-2, join_style=1).buffer(-1e-2, join_style=1)
    # slight rounding. maybe redundant?
    polygon = make_valid(polygon)
    # validate the shape using the shapely built-in method make_valid
    polygon = polygon.difference(Polygon([(-100,-100),(100,-100),(100,0),(-100,0)])) 
    # just to make sure to get rid of the negative region of the polygon
    polygon = unary_union(polygon)
    # re-validate the shape. maybe redundant?

    if (polygon.is_valid & (LineString([(0,0),(1,0)]).intersects(polygon))):
        polygon = translate(polygon, xoff = -polygon.bounds[0])
        polygon = scale(polygon, xfact  = 1.0/(polygon.bounds[2]-polygon.bounds[0])
                               , yfact  = 1.0/(polygon.bounds[3])
                               , origin = (0,0))
        return [polygon, weights_norm]
        # re-scale the intersection-treated DbM shape in a box [0,1] x [0,1]
    else: 
        return [Polygon([(0,0), (1,0), (0,0)]), weights_norm]
        # the root of the rib does not exist, etc. -> return no rib

def getRandomSamplesOnNSphere(nDim, radius, numberofSamples):
    X = np.random.default_rng().normal(size=(numberofSamples, nDim))
    # ref. https://mathworld.wolfram.com/HyperspherePointPicking.html
    X = radius / np.sqrt(np.sum(X**2, 1, keepdims=True)) * X
    return X.tolist()[0]

def lesGridPre(morphedShape, ribWidth, ribHeight, ribSpacing, gridDir):
    global Re_tau
    global N_x
    global N_y
    global N_z
    global L_x
    global L_y
    global L_z

    N_buff = 2 # 2
    N = {'x': N_x, 'y': N_y + 2*N_buff, 'z': N_z}

    L_onePattern = np.round((ribWidth+ribSpacing)/(Re_tau), 15)
    N_r = L_z // L_onePattern
    if 2.*(L_z - N_r * L_onePattern) >= L_onePattern :
        N_r += 1
    L_z = N_r * L_onePattern # rib adjustment

    L_y += (morphedShape.area * (ribWidth/(Re_tau)) * (ribHeight/(Re_tau)) * float(N_r)) / L_z
    # L_y *= (1 + morphedShape.length - morphedShape.exterior.intersection(LineString([(0,0),(1,0)])).length) \
    #         / (2. - morphedShape.area * ribHeight / (Re_tau))

    L_buff = np.round(2./(Re_tau), 15) # np.round(2./(Re_tau), 15)
    L = {'x': L_x, 'y': L_y + 2*L_buff, 'z': L_z}

    y_rib_t = 30. / (Re_tau) # ribHeight/(Re_tau) # top end y coord for a single rib pattern
    n_rib_t = 1 + 32         # 1   + 32           # top end gridpoint for a single rib pattern

    y_hmin = y_rib_t * (Re_tau) / (n_rib_t - 1) * 1.2 # nearest wall y+ in the hyperbolic tangent grid 
    y_hfac = .1; eps = np.inf 
    while (eps > 5e-5): # iterative hyperbolic tangent factor finder for grid sizing
        tmp = (1. - y_hmin/(Re_tau)/((L_y - 2*y_rib_t)*.5))*np.tanh(y_hfac*.5)
        y_hfac_new = np.arctanh(tmp)/(.5-1./(N_y - 2*n_rib_t))
        eps = np.abs(y_hfac_new - y_hfac)
        y_hfac = y_hfac_new

    # begin to set the grid parameters, following the input rule of the LES-IBM code
    gridlines = []
    gridlines.append({'dir'     :                  'X', 'coord_i' :              -L_x/2., 'coord_f'  :        L_x/2., \
                      'index_i' :                    1, 'index_f' :                  N_x, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                   0.})   
           
    gridlines.append({'dir'     :                  'Y', 'coord_i' :              -L_buff, 'coord_f'  :             0, \
                      'index_i' :                    1, 'index_f' :             1+N_buff, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                   0.})       
           
    gridlines.append({'dir'     :                  'Y', 'coord_i' :                   0., 'coord_f'  :       y_rib_t, \
                      'index_i' :             1+N_buff, 'index_f' :       N_buff+n_rib_t, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                   0.})   
    gridlines.append({'dir'     :                  'Y', 'coord_i' :              y_rib_t, 'coord_f'  :   L_y-y_rib_t, \
                      'index_i' :       N_buff+n_rib_t, 'index_f' : N_buff+N_y-n_rib_t+1, 'grid_opt' :           'H', \
                      'factor1' :                  0.5, 'factor2' :               y_hfac})   
    gridlines.append({'dir'     :                  'Y', 'coord_i' :          L_y-y_rib_t, 'coord_f'  :           L_y, \
                      'index_i' : N_buff+N_y-n_rib_t+1, 'index_f' :           N_buff+N_y, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                    0})  
          
    gridlines.append({'dir'     :                  'Y', 'coord_i' :                  L_y, 'coord_f'  :    L_y+L_buff, \
                      'index_i' :           N_buff+N_y, 'index_f' :         2*N_buff+N_y, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                   0.})    
           
    gridlines.append({'dir'     :                  'Z', 'coord_i' :                   0., 'coord_f'  :           L_z, \
                      'index_i' :                    1, 'index_f' :                  N_z, 'grid_opt' :           'U', \
                      'factor1' :                   0., 'factor2' :                   0.})   

    # n_rib_w = 1   + 32 # horizontal gridpoints for a single rib pattern
    # if ribSpacing/ribWidth < .20:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth < .17:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth < .14:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth < .10:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth < .07:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth < .04:
    #     n_rib_w += 1
    # if ribSpacing/ribWidth == 0.:
    #     n_rib_w += 2
    # # gradual adjustment of the # of gridpoints for spacings
    # # if no spacing, then all gridpoints belong to the ribs

    # if ribSpacing != 0.:
    #     z_hmin = ribWidth/(n_rib_w-1) # minimum z+ in the spacing
    #     z_hfac = .1; eps = np.inf 
    #     while (eps > 5e-5): # iterative geometric factor finder for grid sizing
    #         tmp = (1. - z_hmin/(Re_tau)/(((ribSpacing/2)/(Re_tau))*1.))*np.tanh(z_hfac*1.)
    #         z_hfac_new = np.arctanh(tmp)/(1.-1./(((N_z - 1) - (n_rib_w - 1) * N_r) / N_r / 2))
    #         eps = np.abs(z_hfac_new - z_hfac)
    #         z_hfac = z_hfac_new

    # for i in range(N_r): # N_r rib patterns are repeated in the spanwise direction
    #     # left spacing
    #     zcrd_i = np.round(i     * (ribWidth+ribSpacing)/(Re_tau), 15)
    #     zcrd_f = np.round(zcrd_i + (ribSpacing/2)/(Re_tau), 15)
    #     zidx_i = int(1 + i     * ((n_rib_w - 1) + ((N_z - 1) - (n_rib_w - 1) * N_r) / N_r))
    #     zidx_f = int(zidx_i + ((N_z - 1) - (n_rib_w - 1) * N_r) / N_r / 2)
    #     if ribSpacing/ribWidth >= .08:
    #         gridlines.append({'dir'     :     'Z', 'coord_i' :  zcrd_i, 'coord_f'  :  zcrd_f, \
    #                           'index_i' :  zidx_i, 'index_f' :  zidx_f, 'grid_opt' :     'H', \
    #                           'factor1' :       0, 'factor2' : z_hfac})
    #     elif ribSpacing/ribWidth > 0.:
    #         gridlines.append({'dir'     :     'Z', 'coord_i' :  zcrd_i, 'coord_f'  :  zcrd_f, \
    #                           'index_i' :  zidx_i, 'index_f' :  zidx_f, 'grid_opt' :     'U', \
    #                           'factor1' :       0, 'factor2' :     0.})

    #     # rib
    #     zcrd_i = np.round(zcrd_f, 15)
    #     zcrd_f = np.round(zcrd_i + ribWidth/(Re_tau), 15)
    #     zidx_i = int(zidx_f)
    #     zidx_f = int(zidx_i + (n_rib_w - 1))
    #     gridlines.append({'dir'     :     'Z', 'coord_i' :  zcrd_i, 'coord_f'  :  zcrd_f, \
    #                       'index_i' :  zidx_i, 'index_f' :  zidx_f, 'grid_opt' :     'U', \
    #                       'factor1' :       0, 'factor2' :      0.})
    #     # right spacing
    #     zcrd_i = np.round(zcrd_f, 15)
    #     zcrd_f = np.round(zcrd_i + (ribSpacing/2)/(Re_tau), 15) 
    #     zidx_i = int(zidx_f)
    #     zidx_f = int(zidx_i + ((N_z - 1) - (n_rib_w - 1) * N_r) / N_r / 2)
    #     if ribSpacing/ribWidth >= .08:
    #         gridlines.append({'dir'     :     'Z', 'coord_i' :  zcrd_i, 'coord_f'  :  zcrd_f, \
    #                           'index_i' :  zidx_i, 'index_f' :  zidx_f, 'grid_opt' :     'H', \
    #                           'factor1' :       1, 'factor2' :  z_hfac})
    #     elif ribSpacing/ribWidth > 0.:
    #         gridlines.append({'dir'     :     'Z', 'coord_i' :  zcrd_i, 'coord_f'  :  zcrd_f, \
    #                           'index_i' :  zidx_i, 'index_f' :  zidx_f, 'grid_opt' :     'U', \
    #                           'factor1' :       0, 'factor2' :      0.})

    if abs(gridlines[-1]['coord_f'] - L_z) < 5.e-15:
        gridlines[-1]['coord_f'] = L_z 
    # just to avoid machine-precision round-off error ...

    if not(gridValidity(gridlines, N, L)):
        raise RuntimeError('(lesGridPre) input grid information is invalid')
    # validity check

    gridGenSave(gridlines, N, L, gridDir, debug=True)
    # generate gridfile 'grid.out' and save it into gridDir
    # for the sake of debugging. set debug to false in the production stage. it creates a yz-plane meshgrid information in the current folder 

    return None

def gridValidity(gridlines, N, L):
    gridlines_x = list(); gridlines_y = list(); gridlines_z = list()
    i = 0
    for gridline in gridlines:
        i = i+1
        if gridline['dir'] in ['X','x']:
            gridlines_x.append(gridline)
        elif gridline['dir'] in ['Y','y']:
            gridlines_y.append(gridline)
        elif gridline['dir'] in ['Z','z']:
            gridlines_z.append(gridline)
        else:
            print('[Error] Unavailable grid direction at %d. check the direction(X,Y,Z) again.' %(i))
            raise RuntimeError('(gridValidity) grid validity check failed')

    error = list()
    lib_gridvalid.isValid(gridlines_x, N['x'], L['x'], error)
    lib_gridvalid.isValid(gridlines_y, N['y'], L['y'], error)
    lib_gridvalid.isValid(gridlines_z, N['z'], L['z'], error)

    if error:
        for err in error:
            print(err)
        return False

    return True

def gridGenSave(gridlines, N, L, gridDir, debug=False):
    if gridDir[-1] != '/':
        gridDir += '/'

    x_coord = dict(); y_coord = dict(); z_coord = dict() # grid intervals in x,y,z directions 

    # Call the grid functions from the 'gridfunc' library.
    for i, gridline in enumerate(gridlines):
        if gridline['grid_opt'] in ['U','u']: # unifrom gridlines
            lib_gridfunc.uniform(x_coord,y_coord,z_coord,gridline)
        elif gridline['grid_opt'] in ['G','g']: # geometric progression gridlines
            lib_gridfunc.geometric(x_coord,y_coord,z_coord,gridline)
        elif gridline['grid_opt'] in ['H','h']: # hyperbolic tangent gridlines
            lib_gridfunc.hypertan(x_coord,y_coord,z_coord,gridline)
        else:
            print('[Error] Unavailable grid option at line %d. only U, G, H are recognized.' %(i+1))
            raise RuntimeError('(gridGenSave) grid generation interrupted due to input problem')

    with open(gridDir + 'grid.out', 'w') as f:
        f.write('%23d %23d %23d\n' %(N['x'], N['y'], N['z']))
        f.write('%23.15f %22.15f %22.15f\n' %(L['x'], L['y'], L['z']))
        for i in range(1,N['x']+1):
            f.write('%23.15f' %(x_coord[i]))
        f.write('\n')
        for i in range(1,N['y']+1):
            f.write('%23.15f' %(y_coord[i]))
        f.write('\n')
        for i in range(1,N['z']+1):
            f.write('%23.15f' %(z_coord[i]))

    if debug:
        lib_griddebug.plane_grid(y_coord, z_coord, N['y'], N['z'], 'yz', gridDir)

    return None

def ibmBodyPre(morphedShape, ribWidth, ribHeight, ribSpacing, gridDir, ibmFileDir):
    global Re_tau
    xprdic = 'ON'
    yprdic = 'OFF'
    zprdic = 'ON'
    ibmint = 'ON'
    debugopt = {'U_surf_3D' : 'ON',  \
                'V_surf_3D' : 'OFF', \
                'W_surf_3D' : 'OFF'}

    grid = lib_setgrid.setgrid(gridDir, xprdic, yprdic, zprdic)

    with open(ibmFileDir+'ibmpre_prdic.bin', 'w') as f:
        f.write('%d %d %d %d' %(1, 0, 1, 1))

    N_x = grid['N_x']; N_y = grid['N_y']; N_z = grid['N_z']
    Cell_x = grid['Cell_x']; Cell_y = grid['Cell_y']; Cell_z = grid['Cell_z']
    L_x = grid['L_x']; L_y = grid['L_y']; L_z = grid['L_z']

    X = []; Y = []; Z = []
    XM = []; YM = []; ZM = []

    for i in range(0,N_x+1):
        X.append(grid['grid_info_x'][i]['coord'])
        XM.append(grid['grid_info_x'][i]['center'])
    for j in range(0,N_y+1):
        Y.append(grid['grid_info_y'][j]['coord'])
        YM.append(grid['grid_info_y'][j]['center'])
    for k in range(0,N_z+1):
        Z.append(grid['grid_info_z'][k]['coord'])
        ZM.append(grid['grid_info_z'][k]['center'])

    # find inner body pts
    X = np.array(X, order='F'); Y = np.array(Y, order='F'); Z = np.array(Z, order='F');
    XM = np.array(XM, order='F'); YM = np.array(YM, order='F'); ZM = np.array(ZM, order='F');

    nbody = {'u':0, 'v':0, 'w':0} # number of pts in the body defined
    inout = {'u':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F'), \
             'v':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F'), \
             'w':np.empty((N_x+1,N_y+1,N_z+1), dtype=int, order='F')}   # whether pts is in the body or not

    nbody['u'], inout['u'] = find_inout(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,X,YM,ZM,.0)
    nbody['v'], inout['v'] = find_inout(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,XM,Y,ZM,.0)
    nbody['w'], inout['w'] = find_inout(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,XM,YM,Z,.0)

    ninner = {'u':0, 'v':0, 'w':0}                                # Number of inner pts
    fcp = {'u':np.empty([nbody['u'],3], dtype=int, order='F'), \
           'v':np.empty([nbody['v'],3], dtype=int, order='F'), \
           'w':np.empty([nbody['w'],3], dtype=int, order='F')}    # Index of the forcing pts

    nintp = {'u':0, 'v':0, 'w':0}                                 # Number of bdy pts for interpolation
    intptype = {'u':np.empty([nbody['u'],1], dtype=int, order='F'), \
                'v':np.empty([nbody['v'],1], dtype=int, order='F'), \
                'w':np.empty([nbody['w'],1], dtype=int, order='F')}    # Type of the intp. forcing pts
                # (0 : inner, 1 : face, 2 : edge, 3 : single volume)
    intpindx = {'u':np.empty([nbody['u'],3], dtype=int, order='F'), \
                'v':np.empty([nbody['v'],3], dtype=int, order='F'), \
                'w':np.empty([nbody['w'],3], dtype=int, order='F')}    # Intp. forcing pts direction indicator

    geomfac = {'u':np.empty([nbody['u'],3,3,3], order='F'), \
               'v':np.empty([nbody['v'],3,3,3], order='F'), \
               'w':np.empty([nbody['w'],3,3,3], order='F')}   # Geometric factor for interpolation

    # Find pts for interpolation (located at the surface of the body)
    # print('\n*** INTERPOLATION MODE ON ***')
    ufix_x = list(); ufix_y = list(); ufix_z = list()
    lfix_x = list(); lfix_y = list(); lfix_z = list()
    for i in range(1,N_x):
        ufix_x.append(grid['grid_info_x'][i]['upperfix'])
        lfix_x.append(grid['grid_info_x'][i]['lowerfix'])
    for j in range(1,N_y):
        ufix_y.append(grid['grid_info_y'][j]['upperfix'])
        lfix_y.append(grid['grid_info_y'][j]['lowerfix'])
    for k in range(1,N_z):
        ufix_z.append(grid['grid_info_z'][k]['upperfix'])
        lfix_z.append(grid['grid_info_z'][k]['lowerfix'])

    [nintp['u'], ninner['u'], fcp['u'], intptype['u'], intpindx['u']] = \
         lib_ibm_body.findbdy_intp(1, nbody['u'], inout['u'], ufix_x, ufix_y, ufix_z, lfix_x, lfix_y, lfix_z)
    [nintp['v'], ninner['v'], fcp['v'], intptype['v'], intpindx['v']] = \
         lib_ibm_body.findbdy_intp(2, nbody['v'], inout['v'], ufix_x, ufix_y, ufix_z, lfix_x, lfix_y, lfix_z)
    [nintp['w'], ninner['w'], fcp['w'], intptype['w'], intpindx['w']] = \
         lib_ibm_body.findbdy_intp(3, nbody['w'], inout['w'], ufix_x, ufix_y, ufix_z, lfix_x, lfix_y, lfix_z)

    intptype['u'] = intptype['u'][0:nintp['u'],:].copy()    # cut the unused array to nintp
    intptype['v'] = intptype['v'][0:nintp['v'],:].copy()    # (unused) = (nbody) - (nintp)
    intptype['w'] = intptype['w'][0:nintp['w'],:].copy()

    intpindx['u'] = intpindx['u'][0:nintp['u'],:].copy()    # cut the unused array to nintp
    intpindx['v'] = intpindx['v'][0:nintp['v'],:].copy()    # (unused) = (nbody) - (nintp)
    intpindx['w'] = intpindx['w'][0:nintp['w'],:].copy()

    geomfac['u'] = geomfac['u'][0:nintp['u'],:,:,:].copy()  # cut the unused array to nintp
    geomfac['v'] = geomfac['v'][0:nintp['v'],:,:,:].copy()  # (unused) = (nbody) - (nintp)
    geomfac['w'] = geomfac['w'][0:nintp['w'],:,:,:].copy()

    fcp['u'] = fcp['u'][0:nintp['u']+ninner['u'],:].copy()  # cut the unused array to (nintp+ninner)
    fcp['v'] = fcp['v'][0:nintp['v']+ninner['v'],:].copy()
    fcp['w'] = fcp['w'][0:nintp['w']+ninner['w'],:].copy()

    XX = lib_ibm_body.geomfac_preset(X, XM, xprdic)
    YY = lib_ibm_body.geomfac_preset(Y, YM, yprdic)
    ZZ = lib_ibm_body.geomfac_preset(Z, ZM, zprdic)

    geomfac['u'] = lib_ibm_body.geomfac_intp(XX[:,0], YY[:,1], ZZ[:,2], fcp['u'], intpindx['u'], .0)
    geomfac['v'] = lib_ibm_body.geomfac_intp(XX[:,1], YY[:,2], ZZ[:,0], fcp['v'], intpindx['v'], .0)
    geomfac['w'] = lib_ibm_body.geomfac_intp(XX[:,2], YY[:,0], ZZ[:,1], fcp['w'], intpindx['w'], .0)

    # lib_ibmpre.GFIDebug(geomfac, intptype, intpindx, fcp, X, Y, Z)

    nzero = 0
    iszero = np.empty((N_x-1,N_y-1,N_z-1), dtype=int, order='F')

    nzero, iszero = find_zero_nu_sgs(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,XM,YM,ZM,.0)

    ibmPresetSave(xprdic,yprdic,zprdic,ibmint,ibmFileDir)
    ibmDataSave(grid,inout,nintp,ninner,fcp,intpindx,geomfac,nzero,iszero,ibmFileDir,debug=(debugopt['U_surf_3D'] == 'ON'))

    return None

def find_inout(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,X,Y,Z,T):
    global L_y
    nbody = 0
    inout = np.ones((X.size,Y.size,Z.size), dtype=int, order='F')
    z_norm = np.fmod((Z - (ribSpacing/2)/(Re_tau)), (ribSpacing + ribWidth)/(Re_tau)) / (ribWidth/(Re_tau))
    y_norm = Y / (ribHeight/(Re_tau))
    for k in range(Z.size):
        for j in range(int(np.round(Y.size/2))): # range(Y.size)
            # for i in range(X.size):
            if morphedShape.contains(Point((z_norm[k], y_norm[j]))):
                nbody += X.size
                inout[:,j,k] = 0
            elif Y[j] <= 1.E-10 :
                nbody += X.size
                inout[:,j,k] = 0
        for j in range(int(np.round(Y.size/2))+1,Y.size):
            if Y[j] >= L_y - 1.E-10: # 2.+ (morphedShape.area * (ribWidth/(Re_tau)) * (ribHeight/(Re_tau)) * float(N_r)) / np.round((ribWidth+ribSpacing)/(Re_tau) * float(N_r), 15):
                nbody += X.size
                inout[:,j,k] = 0
    if nbody == 0:
        raise RuntimeError('(find_inout) no rib detected - possibly out of valid morphing range')
    return nbody, inout

def find_zero_nu_sgs(morphedShape,ribWidth,ribHeight,ribSpacing,Re_tau,X,Y,Z,T):
    global L_y
    nzero = 0
    iszero = np.ones((X.size-2,Y.size-2,Z.size-2), dtype=int, order='F')
    z_norm = np.fmod((Z - ((ribSpacing/2)/(Re_tau))), (ribSpacing + ribWidth)/(Re_tau)) / (ribWidth/(Re_tau)) 
    y_norm = Y / (ribHeight/(Re_tau))
    for k in range(1,Z.size-2):
        for j in range(1,int(np.round(Y.size/2-1))): # range(1,Y.size-2)
            # for i in range(1,X.size-2):
            if ( morphedShape.contains(Point((z_norm[k  ], y_norm[j  ]))) | \
                 morphedShape.contains(Point((z_norm[k-1], y_norm[j  ]))) | \
                 morphedShape.contains(Point((z_norm[k+1], y_norm[j  ]))) | \
                 morphedShape.contains(Point((z_norm[k  ], y_norm[j-1]))) | \
                 morphedShape.contains(Point((z_norm[k  ], y_norm[j+1]))) ):
                nzero += X.size-2
                iszero[:,j,k] = 0
            elif ( (Y[j]<=1.E-10) | (Y[j-1]<=1.E-10)):
                nzero += X.size-2
                iszero[:,j,k] = 0
        for j in range(int(np.round(Y.size/2-1))+1,Y.size-2):
            if ( (Y[j]>= L_y - 1.E-10) #2. + (morphedShape.area * (ribWidth/(Re_tau)) * (ribHeight/(Re_tau)) * float(N_r)) / np.round((ribWidth+ribSpacing)/(Re_tau) * float(N_r), 15)) \
               | (Y[j+1]>= L_y - 1.E-10) ):#2. + (morphedShape.area * (ribWidth/(Re_tau)) * (ribHeight/(Re_tau)) * float(N_r)) / np.round((ribWidth+ribSpacing)/(Re_tau) * float(N_r), 15))):
                nzero += X.size-2
                iszero[:,j,k] = 0
    return nzero, iszero

def ibmPresetSave(xprdic,yprdic,zprdic,ibmint,ibmFileDir):
    if ibmFileDir[-1] != '/':
        ibmFileDir += '/'

    if xprdic == 'ON':
        xprdic_ = 1
    else:
        xprdic_ = 0
    if yprdic == 'ON':
        yprdic_ = 1
    else:
        yprdic_ = 0
    if zprdic == 'ON':
        zprdic_ = 1
    else:
        zprdic_ = 0
    if ibmint == 'ON':
        ibmint_ = 1
    else:
        ibmint_ = 0

    with open(ibmFileDir+'ibmpre_prdic.bin', 'w') as f:
        f.write('%d %d %d %d' %(xprdic_, yprdic_, zprdic_, ibmint_))

    return None

def ibmDataSave(grid,inout,nintp,ninner,fcp,intpindx,geomfac,nzero,iszero,ibmFileDir,debug=False):
    if ibmFileDir[-1] != '/':
        ibmFileDir += '/'

    lib_ibmpre.grid_preprocessing_data(grid,ibmFileDir)
    lib_ibmpre.ibm_preprocessing_data(nintp,ninner,fcp,intpindx,geomfac,ibmFileDir)
    lib_ibmpre.les_preprocessing_data(nzero,iszero,ibmFileDir)

    if debug==True:
        X = []; Y = []; Z = []
        XM = []; YM = []; ZM = []
        N_x = grid['N_x']; N_y = grid['N_y']; N_z = grid['N_z']
        for i in range(0,N_x+1):
            X.append(grid['grid_info_x'][i]['coord'])
            XM.append(grid['grid_info_x'][i]['center'])
        for j in range(0,N_y+1):
            Y.append(grid['grid_info_y'][j]['coord'])
            YM.append(grid['grid_info_y'][j]['center'])
        for k in range(0,N_z+1):
            Z.append(grid['grid_info_z'][k]['coord'])
            ZM.append(grid['grid_info_z'][k]['center'])
        lib_ibmpre.surf3D(X,YM,ZM,inout['u'], 'u', ibmFileDir)

    return None
