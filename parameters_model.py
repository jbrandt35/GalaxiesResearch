#snapshot parameters

snapshot_num = 125
snapnum_str = '{:04d}'.format(snapshot_num)

filterFile = 'F475W'
hydro_dir = '/storage/home/hhive1/ssethuram6/data/Work/testData/SG64-2020/DD' + snapnum_str + '/'

snapshot_name = 'output_'+snapnum_str

#where the files should go
PD_output_dir = "/storage/home/hhive1/jbrandt35/data/GalaxiesResearch/Practice"
Auto_TF_file = 'snap'+snapnum_str+filterFile + '.logical'
Auto_dustdens_file = 'snap'+snapnum_str+filterFile + '.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/' +snapshot_name+filterFile + '.rtin'
outputfile = PD_output_dir+'/'+snapshot_name+filterFile + '.rtout'


#===============================================
#GRID POSITIONS
#===============================================
x_cent = 0.49448714159030444
y_cent = 0.5108040201005025
z_cent = 0.49365947383978714

TCMB = 2.73
