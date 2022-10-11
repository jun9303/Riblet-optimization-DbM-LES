# Grid input validity checker module
import sys

def isValid(gridlines_i, totalline, totaldomain, error):
    lnum = len(gridlines_i)

    if abs((gridlines_i[lnum - 1]['coord_f'] - gridlines_i[0]['coord_i']) - totaldomain) > 5e-15:
        error.append('[Error] Computational domain length in %s-direction does not match with the current gridline inputs.' \
                     %(gridlines_i[0]['dir'].lower()))
    if gridlines_i[0]['index_i'] != 1:
        error.append('[Error] The first index number is not equal to 1 in %s-direction.' \
                     %(gridlines_i[0]['dir'].lower()))
    if gridlines_i[lnum-1]['index_f'] != totalline:
        error.append('[Error] The final index number does not match with the total gridline number in %s-direction.' \
                     %(gridlines_i[0]['dir'].lower()))

    if lnum > 1:
        for n in range(1,lnum-1):
            if abs(gridlines_i[n]['coord_i'] - gridlines_i[n-1]['coord_f']) > 5e-15:
                error.append('[Error] The inital coordinate value of line %d in %s-direction is not equal to the final coordinate value of the previous line.' \
                             %(n, gridlines_i[0]['dir'].lower()))
            if gridlines_i[n]['index_i'] != gridlines_i[n-1]['index_f']:
                error.append('[Error] The inital index of line %d in %s-direction is not equal to the final index of the previous line.' \
                             %(n, gridlines_i[0]['dir'].lower()))

    if error:
        return False
    else:
        return True

if __name__ == '__main__':
    print('[Error] This is a module file. Please execute run.py')
    sys.exit(1)