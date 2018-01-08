import numpy as np
import astropy.units as u

from synthesizAR.atomic import EmissionModel, Element, list_elements

# Build emission model
temperature = 10.**(np.arange(4.5,8,0.05))*u.K
density = np.logspace(7,11,15)/(u.cm**3)
# Include all ions in CHIANTI
selected_elements = ['calcium','iron','magnesium','nickel','oxygen','silicon','sulfur']
ions = [Element(el, temperature, ion_kwargs={'abundance_filename': 'sun_coronal_1992_feldman'}) 
        for el in selected_elements]
em_model = EmissionModel(density, *ions)

# Calculate emissivity and store in a table
em_model.calculate_emissivity('/storage-home/w/wtb2/data/timelag_synthesis/emissivity_table_dominant.h5',notebook=False)

# Save the whole model
em_model.save('/storage-home/w/wtb2/data/timelag_synthesis/base_emission_model_dominant.json')

