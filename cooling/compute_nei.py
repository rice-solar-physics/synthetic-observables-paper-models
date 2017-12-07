import distributed

import synthesizAR
from synthesizAR.atomic import EmissionModel
from synthesizAR.model_ext import EbtelInterface

client = distributed.Client('tcp://128.42.128.76:8786')

# Restore field
field = synthesizAR.Skeleton.restore('/storage-home/w/wtb2/data/timelag_synthesis/cooling/field_checkpoint/')
# Restore emission model from base
em_model = EmissionModel.restore('/storage-home/w/wtb2/data/timelag_synthesis/base_emission_model.json')
# Set up tasks
tasks = em_model.calculate_ionization_fraction(field,
                                              '/storage-home/w/wtb2/data/timelag_synthesis/cooling/full/ionization_fractions.h5',
                                               interface=EbtelInterface)
# Compute
tasks.compute()
# Save emission model
em_model.save('/storage-home/w/wtb2/data/timelag_synthesis/cooling/full/emission_model.json')

client.close()