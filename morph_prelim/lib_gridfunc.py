# Grid function module
import sys
from numpy import tanh, round

def uniform(x,y,z,line):
    if line['dir'] in ['X','x']: # x-direction gridlines
        xyz = x
    elif line['dir'] in ['Y','y']: # y-direction gridlines
        xyz = y
    elif line['dir'] in ['Z','z']: # z-direction gridlines
        xyz = z
    else:
        print('[Error] Inappropriate direction option is implemented.')
        sys.exit(1)

    delta = (line['coord_f'] - line['coord_i']) / (line['index_f'] - line['index_i'])

    for i in range(line['index_i'], line['index_f']+1):
        xyz[i] = line['coord_i']+(i-line['index_i'])*delta

    print('%s-dir., Uniform,   %15.11f ~ %15.11f, Interval is %.6f.'
           %(line['dir'].lower(), line['coord_i'], line['coord_f'], delta))

def geometric(x,y,z,line):
    if line['dir'] in ['X','x']: # x-direction gridlines
        xyz = x
    elif line['dir'] in ['Y','y']: # y-direction gridlines
        xyz = y
    elif line['dir'] in ['Z','z']: # z-direction gridlines
        xyz = z
    else:
        print('[Error] Inappropriate direction option is implemented.')
        sys.exit(1)

    factor1 = line['factor1'] # Expansion(1) or Compression(0)
    factor2 = line['factor2'] # Geometric Ratio (> 1)

    if factor2 <= 1:
        print('[Error] For geometric progression, ratio must be larger than 1')
        sys.exit(1)
    if factor1 != 1:
        factor2 = 1/factor2

    initdelta = (line['coord_f']-line['coord_i'])*(1-factor2)/(1-factor2**(line['index_f']-line['index_i']))
    finaldelta = initdelta*factor2**(line['index_f']-line['index_i']-1)

    for i in range(line['index_i'], line['index_f']+1):
        xyz[i] = line['coord_i']+initdelta*(1-factor2**(i-line['index_i']))/(1-factor2)

    print('%s-dir., Geometric, %15.11f ~ %15.11f, Interval from %.6f to %.6f.'
           %(line['dir'].lower(), line['coord_i'], line['coord_f'], initdelta, finaldelta))


def hypertan(x,y,z,line):
    if line['dir'] in ['X','x']: # x-direction gridlines
        xyz = x
    elif line['dir'] in ['Y','y']: # y-direction gridlines
        xyz = y
    elif line['dir'] in ['Z','z']: # z-direction gridlines
        xyz = z
    else:
        print('[Error] Inappropriate direction option is implemented.')
        sys.exit(1)

    factor1 = line['factor1'] # Expansion(1) or Compression(0) or Symmetric(0.5)
    factor2 = line['factor2'] # Gamma value

    if factor1 != 0:
        for i in range(line['index_i'], line['index_f']+1):
            prop = (i-line['index_i'])/(line['index_f']-line['index_i'])
            xyz[i] = line['coord_i']+(line['coord_f']-line['coord_i'])*factor1* \
                     (1-tanh(factor2*(factor1-prop))/tanh(factor2*factor1))
    else:
        factor1 = 1
        for i in range(line['index_i'], line['index_f']+1):
            prop = (i-line['index_i'])/(line['index_f']-line['index_i'])
            xyz[i] = line['coord_i']+(line['coord_f']-line['coord_i'])*factor1* \
                     tanh(factor2*prop)/tanh(factor2*factor1)


    initdelta = xyz[line['index_i']+1] - xyz[line['index_i']]
    middelta = xyz[round((line['index_i']+line['index_f'])/2)+1] - xyz[round((line['index_i']+line['index_f'])/2)]
    finaldelta = xyz[line['index_f']] - xyz[line['index_f']-1]

    print('%s-dir., Hypertan,  %15.11f ~ %15.11f, Interval from %.6f through %.6f to %.6f.'
           %(line['dir'].lower(), line['coord_i'], line['coord_f'], initdelta, middelta, finaldelta))

if __name__ == '__main__':
    print('[Error] This is a module file. Please execute run.py')
    sys.exit(1)