import astropy.units as u
import distributed

import synthesizAR
from synthesizAR.instruments import InstrumentSDOAIA
from synthesizAR.atomic import EmissionModel

# Connect to client
client = distributed.Client('tcp://127.0.0.1:34705')
# Load field
field = synthesizAR.Skeleton.restore('/storage-home/w/wtb2/data/timelag_synthesis/high_frequency/field_checkpoint/')
# Load emission model
em_model = EmissionModel.restore('/storage-home/w/wtb2/data/timelag_synthesis/high_frequency/nei/emission_model_dominant.json')
# Build instrument and observer
aia = InstrumentSDOAIA([0, 30000]*u.s, use_temperature_response_functions=False)
observer = synthesizAR.Observer(field, [aia], parallel=True)
observer.build_detector_files('/storage-home/w/wtb2/data/timelag_synthesis/high_frequency/nei/',
                              ds=field._convert_angle_to_length(1.0*u.arcsec))
# Build task graph
tasks = observer.flatten_detector_counts(emission_model=em_model)
# Compute tasks
tasks['SDO_AIA_parameters'].compute()
tasks['SDO_AIA_counts'].compute()