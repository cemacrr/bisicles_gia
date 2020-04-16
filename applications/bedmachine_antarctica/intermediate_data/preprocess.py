

EXTERNAL_DATA_PATH='../external_data'
BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineAntarctica_2019-11-05_v01.nc')
MEASURES_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'antarctica_ice_velocity_450m_v2.nc')
MORLIGHEM_TEMPERATURE_NC ='{}/{}'.format(EXTERNAL_DATA_PATH,'AntarcticTemperature-2020-02-14.nc')
INTERMEDIATE_DATA_PATH='.'

def name(s):
    return '{}/antarctica_bedmachine_{}.nc'.format(INTERMEDIATE_DATA_PATH,s)

def temp_name(s,n_layer):
    return '{}/antarctica_bedmachine_temperature_morlighem_{}_{}.nc'.format(
        INTERMEDIATE_DATA_PATH,s,n_layer)

OUTPUT_NC = name('500m')
N_LAYER=24
OUTPUT_TEMP_MOR_NC = temp_name('4km',N_LAYER)

import os, time

def create_new_only(file_name, fun):
    if os.path.isfile(file_name):
        mb = int(os.path.getsize(file_name)/1024/1024)
        mt = time.ctime(os.path.getsize(file_name))
        print ('{} ({} Mb {}) exists - delete if you want to recreate it'.format(file_name, mb, mt ))
    else:
        print ('creating {} ...'.format(file_name))
        fun(file_name)
        print ('...done')

def preprocess_500m(output_nc):        
    from preprocess_thk_bed_btrc import preprocess
    preprocess(output_nc, BEDMACHINE_NC, MEASURES_NC)
    return None
                        
create_new_only(OUTPUT_NC, preprocess_500m)

#create coarse grid versions for later convenience
suffix = ['500m','1km','2km','4km','8km']

for i in range(1,len(suffix)):
    fine = name(suffix[i-1])
    coarse = name(suffix[i])
    def f(file_name):
        from coarsen import coarsenc
        coarsenc(file_name , fine)
        return(None)
    create_new_only(coarse, f)


def preprocess_T_4km(file_name):
    from preprocess_temp import morlighem_temperature_nc, sigma
    sigma_face, sigma_mid = sigma(N_LAYER)
    morlighem_temperature_nc( file_name, name('4km'), MORLIGHEM_TEMPERATURE_NC, sigma_mid)

create_new_only(OUTPUT_TEMP_MOR_NC, preprocess_T_4km)

