#from yt.fields.api import ValidateParameter
from math import pi
from yt.utilities.physical_constants import planck_constant_cgs
import numpy as np

def _metallicity3(field, data):
    return data[('enzo', 'SN_Colour')] / data[('enzo', 'Density')]
add_field(function=_metallicity3,
          name=('gas', 'metallicity3'), units='Zsun',
          display_name='Pop III Metallicity',
          validators=[ValidateDataField(('enzo', 'SN_Colour'))], sampling_type='cell')

def _total_metallicity(field, data):
    if ('enzo', 'SN_Colour') in data.ds.field_list:
        return data[('gas', 'metallicity3')] + data[('gas', 'metallicity')]
    else:
        return data[('gas', 'metallicity')]
add_field(function=_total_metallicity,
          name=('gas', 'total_metallicity'), units='Zsun',
          display_name='Total Metallicity', sampling_type='cell')

def _electron_fraction(field, data):
    return data[('enzo', 'Electron_Density')] / data[('enzo', 'Density')]
add_field(function=_electron_fraction,
          name=('gas', 'electron_fraction'), units='dimensionless',
          display_name='Electron Fraction',
          validators=[ValidateDataField(('enzo', 'Electron_Density'))], sampling_type='cell')

sigma_H2I = YTQuantity(3.71e-18, 'cm**2')
sigma_HI = YTQuantity(6.3e-18, 'cm**2')
J21_norm = YTQuantity(1e21, 'cm**2/erg')
E_LW = YTQuantity(12.4, 'eV')
E_HI = YTQuantity(13.6, 'eV')
nu_H = YTQuantity(3.2881e15, 'Hz')
mh = YTQuantity(1.673e-24, 'g')
def _J21_LW(field, data):
    return J21_norm * data[('enzo', 'H2I_kdiss')] * E_LW / \
        (4.0 * pi*pi * sigma_H2I * nu_H)
add_field(function=_J21_LW, name=('gas', 'J21_LW'),
          units='dimensionless', display_name='J$_{21}$',
          validators=[ValidateDataField(('enzo', 'H2I_kdiss'))], sampling_type='cell')

def _J_LW(field, data):
    return data[('enzo', 'H2I_kdiss')] * E_LW / \
        (4.0 * pi*pi * sigma_H2I * nu_H)
add_field(function=_J_LW,
          name=('gas', 'J_LW'),
          units='erg/cm**2', display_name='J$_{LW}$',
          validators=[ValidateDataField(('enzo', 'H2I_kdiss'))], sampling_type='cell')

def _J_Lyman(field, data):
    return data[('enzo', 'HI_kph')] * planck_constant_cgs / \
        (4.0 * pi*pi * sigma_HI)
add_field(function=_J_Lyman,
          name=('gas', 'J_Lyman'), units='erg/cm**2', display_name='J(Lyman)',
          validators=[ValidateDataField(('enzo', 'HI_kph'))], sampling_type='cell')

def _column_density(field, data):
    nH2 = (-1.0/sigma_HI) * np.log((sigma_H2I/sigma_HI) *
                                   (np.maximum(data[('enzo', 'HI_kph')], YTQuantity(1e-199, '1/s')) /
                                    data[('enzo', 'H2I_kdiss')]))
    return np.maximum(YTQuantity(1, '1/cm**2'), nH2)
add_field(function=_column_density,
          name=('gas', 'column_density'), units='1/cm**2', display_name='Column Density',
          validators=[ValidateDataField(('enzo', 'HI_kph')),
                      ValidateDataField(('enzo', 'H2I_kdiss'))], sampling_type='cell')

def _gas_fraction(field, data):
    return data[('gas', 'density')] / data[('gas', 'matter_density')]
add_field(function=_gas_fraction,
          name=('gas', 'gas_fraction'), units='dimensionless', display_name='Gas Fraction',
          validators=[ValidateDataField(('enzo', 'Dark_Matter_Density'))], sampling_type='cell')

def _dark_matter_mass(field, data):
    return data[('enzo', 'Dark_Matter_Density')] * data[('gas', 'cell_volume')]
add_field(function=_dark_matter_mass, name=('gas', 'dark_matter_mass'), units='g', display_name='Dark Matter Mass',
          validators=[ValidateDataField(('enzo', 'Dark_Matter_Density'))], sampling_type='cell')

def _mass_divergence(field, data):
    return data[('gas', 'cell_mass')] * data[('gas', 'velocity_divergence')]
add_field(function=_mass_divergence, name=('gas', 'mass_divergence'), units='g/s',
          display_name='Mass Divergence',
          validators=[ValidateParameter('center'),
                      ValidateParameter('bulk_velocity')], sampling_type='cell')

def _radial_mass_flux(field, data):
    return -4*np.pi * data[('index', 'radius')]**2 * data[('gas', 'density')] * data[('gas', 'radial_velocity')]
add_field(function=_radial_mass_flux,
          name=('gas', 'radial_mass_flux'), units='g/s',
          display_name='Radial Mass Flux',
          validators=[ValidateParameter('center'),
                      ValidateParameter('bulk_velocity')], sampling_type='cell')

def _inertia_tensor_xx(field, data):
    center = data.get_field_parameter('center')
    dely = data[('index', 'y')] - center[1]
    delz = data[('index', 'z')] - center[2]
    return data[('gas', 'density')] * data[('gas', 'cell_volume')] * (dely**2 + delz**2)
def _inertia_tensor_yy(field, data):
    center = data.get_field_parameter('center')
    delx = data[('index', 'x')] - center[0]
    delz = data[('index', 'z')] - center[2]
    return data[('gas', 'density')] * data[('gas', 'cell_volume')] * (delx**2 + delz**2)
def _inertia_tensor_zz(field, data):
    center = data.get_field_parameter('center')
    delx = data[('index', 'x')] - center[0]
    dely = data[('index', 'y')] - center[1]
    return data['density'] * data[('gas', 'cell_volume')] * (delx**2 + dely**2)
def _inertia_tensor_xy(field, data):
    center = data.get_field_parameter('center')
    delx = data[('index', 'x')] - center[0]
    dely = data[('index', 'y')] - center[1]
    return data[('gas', 'density')] * data[('gas', 'cell_volume')] * delx * dely
def _inertia_tensor_xz(field, data):
    center = data.get_field_parameter('center')
    delx = data[('index', 'x')] - center[0]
    delz = data[('index', 'z')] - center[2]
    return data[('gas', 'density')] * data[('gas', 'cell_volume')] * delx * delz
def _inertia_tensor_yz(field, data):
    center = data.get_field_parameter('center')
    dely = data[('index', 'y')] - center[1]
    delz = data[('index', 'z')] - center[2]
    return data[('gas', 'density')] * data[('gas', 'cell_volume')] * dely * delz

for dim in "xyz":
    add_field(name=('gas', "inertia_tensor_%c%c" % (dim, dim)),
              sampling_type="cell",
              function=eval("_inertia_tensor_%c%c" % (dim, dim)),
              units="g*cm**2",
              validators=[ValidateParameter('center')])
add_field(name=('gas', "inertia_tensor_xy"),
          sampling_type="cell",
          function=_inertia_tensor_xy,
          units="g*cm**2",
          validators=[ValidateParameter('center')])
add_field(name=('gas', "inertia_tensor_xz"),
          sampling_type="cell",
          function=_inertia_tensor_xz,
          units="g*cm**2",
          validators=[ValidateParameter('center')])
add_field(name=('gas', "inertia_tensor_yz"),
          sampling_type="cell",
          function=_inertia_tensor_yz,
          units="g*cm**2",
          validators=[ValidateParameter('center')])

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def _H_nuclei_density(field, data):
    if ("gas", "H_p0_number_density") not in data.ds.field_info:
        return 0.76 * data[("gas", "density")] / mh
    field_data = np.zeros_like(data[("gas", "H_p0_number_density")])
    for species in ["H_p0", "H_p1", "H_m1", "H2_p0", "H2_p1", "HD_p0"]:
        if ("gas", "%s_number_density" % species) not in data.ds.field_info:
            continue
        nucleus = species[:species.find("_")]
        num = ""
        for i in nucleus[nucleus.find("H"):]:
            if i.isdigit():
                num += i
            else:
                break
        if not num:
            num = "1"
        num = int(num)
        field_data += num * data[("gas", "%s_number_density" % species)]
    return field_data
add_field(function=_H_nuclei_density,
          name=("gas", "H_nuclei_density"), units="cm**-3", sampling_type='cell')

def _log_N(field, data):
    return np.log10(data[("gas", "column_density")])
add_field(function=_log_N, name=("gas", "log_N"), sampling_type='cell')

def _log_J(field, data):
    return np.log10(np.maximum(data[("gas", "J_Lyman")], YTQuantity(1e-99, 'erg/cm**2')))
add_field(function=_log_J, name=("gas", "log_J"), sampling_type='cell')

def _log_Z(field, data):
    return np.log10(data[("gas", "metallicity")])
add_field(function=_log_Z, name=("gas", "log_Z"), sampling_type='cell')

def _log_T(field, data):
    return np.log10(data[("gas", "temperature")])
add_field(function=_log_T, name=("gas", "log_T"), sampling_type='cell')

def _log_nH(field, data):
    return np.log10(data[("gas", "H_nuclei_density")])
add_field(function=_log_nH, name=("gas", "log_nH"), sampling_type='cell')


########################################################################
########################################################################
########################################################################

def p3_sn(pfilter, data):
    return (data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') < 1e-10)
add_particle_filter(function=p3_sn, name='p3_sn', requires=['particle_mass', 'particle_type'])

def p3_bh(pfilter, data):
    return (data[('nbody', 'particle_type')] == 1) & (data[('nbody', 'creation_time')] > 0) & \
        (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3)
add_particle_filter(function=p3_bh, name='p3_bh', requires=['creation_time', 'particle_mass', 'particle_type'])

def p3_living(pfilter, data):
    return (data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3)
add_particle_filter(function=p3_living, name='p3_living', requires=['particle_mass', 'particle_type'])

def p3(pfilter, data):
    return ((data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') < 1e-10)) \
        |  ((data[('nbody', 'particle_type')] == 1) & (data[('nbody', 'creation_time')] > 0) & \
            (data[('nbody', 'particle_mass')].in_units('Msun') > 1)) \
        |  ((data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3))
add_particle_filter(function=p3, name='p3',
                requires=['particle_mass', 'particle_type', 'creation_time'])

def p2(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7)
add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])

def p2_young(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7) & (data[('nbody', 'age')].in_units('Myr') <= 20)
add_particle_filter(function=p2_young, name='p2_young', requires=['particle_mass', 'particle_type', 'creation_time'])

def p2_old(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7) & (data[('nbody', 'age')].in_units('Myr') > 20)
add_particle_filter(function=p2_old, name='p2_old', requires=['particle_mass', 'particle_type', 'creation_time'])

def stars(pfilter, data):
    return ((data[('nbody', 'particle_type')] == 2) & (data[('nbody', 'creation_time')] > 0))
add_particle_filter(function=stars, name='stars', requires=['particle_mass', 'particle_type', 'creation_time'])

def dm(pfilter, data):
    return ((data[('nbody', 'particle_type')] == 1) & (data[('nbody', 'creation_time')] < 0))
add_particle_filter(function=dm, name='dm', requires=['particle_type', 'creation_time'])

def finest(pfilter, data):
    return data[('nbody', 'particle_type')] == 4
add_particle_filter(function=finest, name="mrp", requires=['particle_type'])
