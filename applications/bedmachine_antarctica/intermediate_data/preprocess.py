

EXTERNAL_DATA_PATH='../external_data'
BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineAntarctica_2019-11-05_v01.nc')
MEASURES_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'antarctica_ice_velocity_450m_v2.nc')

INTERMEDIATE_DATA_PATH='.'
OUTPUT_NC = '{}/{}'.format(INTERMEDIATE_DATA_PATH,'antarctica_bedmachine_500m.nc')

import os

if os.path.isfile(OUTPUT_NC):
    print ('{} exists - delete if you want to recreate it'.format(OUTPUT_NC))
else:
    from preprocess_thk_bed_btrc import preprocess
    preprocess(OUTPUT_NC, BEDMACHINE_NC, MEASURES_NC)
