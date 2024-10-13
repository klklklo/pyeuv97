import pandas as pd
import xarray as xr
import numpy as np
import pyeuv97._misc as _m


class EUV97:
    def __init__(self):
        self._bands_dataset, self._lines_dataset = _m.get_euv97_coeffs()
        self._bands_coeffs = np.array(
            self._bands_dataset[['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9']].to_dataarray()).T
        self._lines_coeffs = np.array(
            self._lines_dataset[['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9']].to_dataarray()).T

    def prepare_data(self, proxy):
        if isinstance(proxy, tuple):
            return np.hstack([np.append(1, proxy), np.append(1, proxy)])[None, :]
        array = np.array([np.hstack([np.append(1, i), np.append(1, i)]) for i in proxy])
        return array

    def get_spectral_lines(self, proxy):
        x = self.prepare_data(proxy)
        res = np.dot(self._lines_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('line', 'lac'), res)},
                          coords={'line': self._lines_dataset['line'].values,
                                  'lyman_alpha': x[:, 1],
                                  'HeI': x[:, 2],
                                  'F10.7_daily': x[:, 3],
                                  'F10.7_81avg': x[:, 4],
                                  })

    def get_spectral_bands(self, proxy):
        x = self.prepare_data(proxy)

        res = np.dot(self._bands_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('band_center', 'lac'), res),
                                     'lband': ('band_number', self._bands_dataset['lband'].values),
                                     'uband': ('band_number', self._bands_dataset['uband'].values),
                                     'center': ('band_number', self._bands_dataset['center'].values),
                                     'band_width': ('band_number', self._bands_dataset['width'].values)
                                     },
                          coords={'line': self._bands_dataset['line'].values,
                                  'lyman_alpha': x[:, 1],
                                  'HeI': x[:, 2],
                                  'F10.7_daily': x[:, 3],
                                  'F10.7_81avg': x[:, 4],
                                  'band_number' : np.arange(23)
                                  })

    def get_spectra(self, proxy):
        return self.get_spectral_bands(proxy), self.get_spectral_lines(proxy)

e = EUV97()
mas = [(2,3,4,5),(4,5,6,7)]
mas1 = (2,3,4,5)
# print(e.get_spectral_bands(mas1))
print(e.get_spectral_lines(mas))

# print(0 + 0*2+0*3+0*4+0*5 - 1.06994e9 + 0*2+0*3+7.94298e6*4+9.47885e6*5)
# print(-1.16091e9 - -9.34911e-5*2+1.01294e-2*3+0*4+0*5 + 0 + 0*2+0*3+0*4+0*5)
# print(-9.90773830e+08)

# p = pd.read_csv('_coeffs/coeffs1.csv')
# d = p.loc[p['k'] == 2]
# d.to_csv('_coeffs/k2_coeffs.csv', index=False)

# p = pd.read_csv('_coeffs/prepared_data.csv')
# k2 = p[p['k'] == 2][['a0', 'a1', 'a2', 'a3', 'a4']]
# k2.columns = ['a5', 'a6', 'a7', 'a8', 'a9']
# k2 = k2.reset_index(drop=True)
#
# k = p.drop(p[p['k'] == 2].index)
# k = k.reset_index(drop=True)
#
# k_all = pd.concat([k, k2], axis=1)
# k_all = k_all.drop(columns=['k'])
# k_all.to_csv('prepared_intervals.csv', index = False)
# xr.Dataset.from_dataframe(k_all).to_netcdf('prepared_intervals.nc')
# print(xr.open_dataset('prepared_intervals.nc').dtypes)

# k1 = p[p['k'] == 1][['a0', 'a1', 'a2', 'a3', 'a4']]
# k2 = p[p['k'] == 2][['a0', 'a1', 'a2', 'a3', 'a4']]
# k2.columns = ['a5', 'a6', 'a7', 'a8', 'a9']
# k1 = k1.reset_index(drop=True)
# k2 = k2.reset_index(drop=True)
# k_all = pd.concat([k1, k2], axis=1)
# k_all.to_csv('_coeffs/prepared_data.csv', index=False)

# p = pd.read_csv('_coeffs/prepared_data.csv')
# print(p)
# xr.Dataset(p).to_netcdf('_coeffs/prepared_data.nc')
# print(xr.open_dataset('_coeffs/prepared_data.nc'))

# p = pd.read_csv('_coeffs/individ_lines.csv')
# k2 = p[p['k'] == 2][['a0', 'a1', 'a2', 'a3', 'a4']]
# k2.columns = ['a5', 'a6', 'a7', 'a8', 'a9']
# k2 = k2.reset_index(drop=True)
#
# k = p.drop(p[p['k'] == 2].index)
# k = k.reset_index(drop=True)
# k_all = pd.concat([k, k2], axis=1)
# k_all = k_all.drop(columns=['k', 'uband'])
# k_all = k_all.rename(columns={'lband':'line'})
# k_all.to_csv('_coeffs/prepared_lines.csv', index=False)
# xr.Dataset().from_dataframe(k_all).to_netcdf('_coeffs/prepared_lines.nc')