import xesmf
import xarray as xr
from numpy import arange, append

ds = xr.open_dataset('data_byhour_plev.nc')

# Convert from g[OH]/g[air] to molec[OH] / cm3:
# V_air (cm3) = nRT/P * 1.e-6       # V{air]: volume of air (m3); n: moles of air; T: temperature; P: pressure, R: gas constant
# n = g[air] / M[air]               # M[air]: molar mass of air (g/mol)
# N[OH] = g[OH] * Na / M[OH]        # Na = Avogadro Number; M[OH]: molar mass of OH (g/mol); N[OH]: molec of OH
# 
# Conversions:
# 0: initial data in g/g        : x = m_OH / m_air       # g(oh) / g(air)
# 1. convert to mol ratio       : x = x * M_air / M_oh   # mol(oh) / mol(air)
# 2. convert OH to molecules    : x = x * Na             # molec(oh) / mol(air)
# 3. get volume from air moles  : x = x * P / RT         # molec(oh) / m3
# 4. convert to volume in cm3   : x = x * 1.e-6          # molec(oh) / cm3


M_air = 28.964
M_oh = 15.9994 + 1.0079
Na = 6.02205e23
R = 8.3144

# Do the spatial regridding:
lats = arange(-89.5, 90)
lons = arange(-179.5, 180)
glb1x1 = xr.Dataset(coords=dict(lon=lons, lat=lats))
regridder = xesmf.Regridder(ds, glb1x1, method='bilinear')
oh = regridder(ds)

# Calculate OH in the right unit
OH = oh['oh'] * M_air / M_oh * Na * (oh['pressure_level'][:] * 100) / (R * oh['t'][:]) * 1.e-6

# Create fake hybrid coefficients: 1st level follows the surface pressure (a = 0, b = 1), all the others are at set pressures (b = 0)
p = ds['pressure_level'][:].values.astype('float32')
za = (p[1:] + p[:-1]) / 2.
za = append(za, 0)
za = append(0, za)
#za = za[:-1] - za[1:]
zb = za * 0
zb[0] = 1.

# a should be in Pa:
za *= 100

# Store:
oh['z_a_grid'] = ('nlev', za)
oh['z_b_grid'] = ('nlev', zb)

OH.attrs['units'] = 'molec / cm ** 3'
for k, v in ds.oh.attrs.items():
    oh.oh.attrs[k] = v

# Compute monthly averaging
oh['OH'] = OH.groupby(OH.valid_time.dt.month).mean().astype('double').fillna(0)
print(oh.OH.values.shape)

# There's some missing data at lon index 179. Fill it in by interpolation:
oh.OH.values[:, :, :, 179] = (oh.OH.values[:, :, :, 178] + oh.OH.values[:, :, :, 180])/2.

#oh['OH'].values = oh.OH.values[:, ::-1, :, :]
#oh.to_netcdf('oh_cams_2021_v5.nc')
oh.to_netcdf('oh_cams_2021_monthly.nc')
