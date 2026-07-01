import xarray as xr
import sys


def trim_halos_from_emission_file(filename: str) -> None:
    print(filename)
    ds = xr.open_dataset(filename)
    ds.load()
    ds.close()
    df = ds.to_dataframe()

    # Halos cover to one grid cell of the parent grid:
    # glb6x4 -> eur3x2: 2 grid-cells on each side (6/3 and 4/2)
    # eur3x2 -> gns1x1: 3 grid-cells on the lons (3/1) and 2 in lat (2/1)

    df.loc[:, 'is_halo'] = True
    df.loc[df.region == 'glb600x400', 'is_halo'] = False

    # Remove the halos from the eur3x2 region:
    lons_eur = df[df.region == 'eur300x200'].lon.drop_duplicates().values[2:-2]
    lats_eur = df[df.region == 'eur300x200'].lat.drop_duplicates().values[2:-2]
    df.loc[(df.region == 'eur300x200') & df.lon.isin(lons_eur) & df.lat.isin(lats_eur), 'is_halo'] = False

    # Same for the gns1x1 region:
    lons_gns = df[df.region == 'gns100x100'].lon.drop_duplicates().values[3:-3]
    lats_gns = df[df.region == 'gns100x100'].lat.drop_duplicates().values[2:-2]
    df.loc[(df.region == 'gns100x100') & df.lon.isin(lons_gns) & df.lat.isin(lats_gns), 'is_halo'] = False

    ds_out = df[~df.is_halo].drop(columns='is_halo').to_xarray()

    # Re-add the arguments from the original dataset:
    for varname in ds.data_vars:
        for k, v in ds[varname].attrs.items():
            ds_out[varname].attrs[k] = v

    ds_out = ds_out.drop_vars(['ntc', 'ng', 'nmon'])

    for k, v in ds.attrs.items():
        ds_out.attrs[k] = v

    ds_out.to_netcdf(filename)


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        trim_halos_from_emission_file(fname)
