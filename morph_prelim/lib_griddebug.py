def midpoints(x_coord, y_coord, z_coord, N_x, N_y, N_z, debug_dir):
	if debug_dir[-1] != '/':
		debug_dir += '/'
	file = open(debug_dir + 'grid_midp.dat', 'w')
	file.write('XMP\n')
	for i in range(1,N_x):
	    file.write('%23.15f' %((x_coord[i]+x_coord[i+1])/2))
	file.write('\nYMP\n')
	for i in range(1,N_y):
	    file.write('%23.15f' %((y_coord[i]+y_coord[i+1])/2))
	file.write('\nZMP\n')
	for i in range(1,N_z):
	    file.write('%23.15f' %((z_coord[i]+z_coord[i+1])/2))
	file.close()

def plane_grid(i_coord, j_coord, N_i, N_j, planename, debug_dir):
	if debug_dir[-1] != '/':
		debug_dir += '/'
	file = open(debug_dir+('plane_%s.dat' % planename), 'w')
	file.write('VARIABLES="%s","%s"\n' %(planename[0], planename[1]))
	file.write('ZONE I=%d,J=%d,F=POINT\n' %(N_i, N_j))
	for j in range(1,N_j+1):
		for i in range(1,N_i+1):
			file.write('%13.5f %12.5f' %(i_coord[i], j_coord[j]))
			if not ((i==N_i+1) and (j==N_j+1)):
				file.write('\n')
	file.close()

def deltaplot(coord, N, linename, debug_dir):
	if debug_dir[-1] != '/':
		debug_dir += '/'
	file = open(debug_dir + ('del_%s_plot.dat' % linename), 'w')
	file.write('%-5s,%-13s,%-13s\n' %('i','x','dx'))
	for i in range(1,N):
		file.write('%5d,%13.5f,%13.5f' %(i, coord[i], coord[i+1]-coord[i]))
		if not i==N:
			file.write('\n')
	file.close()