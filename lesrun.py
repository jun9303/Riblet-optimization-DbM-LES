import sys, os, stat, subprocess, shutil
import numpy as np
import matplotlib.pyplot as plt

def solverRun(gridDir, ibmDir, testDir):
    if ibmDir[-1] != '/':
        ibmDir += '/'
    if testDir[-1] != '/':
        testDir += '/'

    solverSrc = './solver/solver_exec'
    solverDst = testDir+'solver'
    shutil.copy(solverSrc, solverDst)
    os.chmod(solverDst, stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR| \
                        stat.S_IRGRP|stat.S_IWGRP|stat.S_IXGRP| \
                        stat.S_IROTH|stat.S_IWOTH|stat.S_IXOTH )
    # grant 'chmod 777' permission to a solver executable copy to make its run certain

    settingSrc = './solver/settings.in'
    settingDst = testDir+'settings.in'
    shutil.copy(settingSrc, settingDst)
    os.chmod(settingDst, stat.S_IRUSR|stat.S_IWUSR| \
                         stat.S_IRGRP|stat.S_IWGRP| \
                         stat.S_IROTH|stat.S_IWOTH )
    # grant 'chmod 666' permission to a setting file copy to make its read/write certain

    bcSrc = './solver/boundary.in'
    bcDst = testDir+'boundary.in'
    shutil.copy(bcSrc, bcDst)
    os.chmod(bcDst, stat.S_IRUSR|stat.S_IWUSR| \
                    stat.S_IRGRP|stat.S_IWGRP| \
                    stat.S_IROTH|stat.S_IWOTH )
    # grant 'chmod 666' permission to a setting file copy to make its read/write certain

    # prefldSrc = './solver/fld000000'
    # prefldDst = testDir+'fld000000'
    # shutil.copy(prefldSrc, prefldDst)
    # os.chmod(prefldDst, stat.S_IRUSR|stat.S_IWUSR| \
    #                     stat.S_IRGRP|stat.S_IWGRP| \
    #                     stat.S_IROTH|stat.S_IWOTH )
    # # grant 'chmod 666' permission to a prefld file copy to make its read/write certain

    cwd_ori = os.getcwd()

    os.chdir(testDir)

    os.mkdir('output')
    os.mkdir('output/ftr')
    os.mkdir('output/field')
    os.mkdir('output/field_avg')
    # create output folders

    try:
        subprocess.run(['./solver', gridDir, ibmDir], check=True, text=True)
    except subprocess.CalledProcessError as e:
        print('(solverRun) subprocess error: \n' + str(e))
        return -0.1

    # try:
    #     result = float(result.stdout)
    # except:
    #     print('(solverRun) not a numeric output in stdout')
    #     return -1.0

    os.chdir(cwd_ori)

    return None
