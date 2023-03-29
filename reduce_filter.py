import numpy as np

number_of_wavelengths_requested = 10

filter = np.loadtxt("JWST_test.filter")

new_filter = filter[::int(len(filter)/number_of_wavelengths_requested)]

np.savetxt(f"JWST_test_{number_of_wavelengths_requested}.filter", new_filter, delimiter = " ")


