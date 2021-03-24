# Determinant-based-Fast-Greedy-Sensor-Selection-Algorithm
This repository contains Matlab (R2013a) code to reproduce results for the Determinant-based Fast Greedy Sensor Selection Algorithm.

Due to GitHub file size limitations, a dataset is linked online:[NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)
- sst.wkmean.1990-present.nc
- lsmask.nc


**Ocean sea surface temperature data is provided by NOAA.
NOAA_OI_SST_V2 data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at https://www.esrl.noaa.gov/psd/.**

## License
[MIT-License](https://github.com/YujiSaitoJapan/Determinant-based-Fast-Greedy-Sensor-Selection-Algorithm/blob/add-license-1/LICENSE)

## Code
- main.m
  - read_NOAA_SST.m (reading data)
  - make_ensopod.m (truncating data)
 
 (sensor selection)
  - sensor_random.m
  - sensor_convex.m
    - subsensor_convex.m
      - sens_sel_approxnt_vec
  - sensor_QR.m
    - subsensor_QR.m
  - sensor_DG.m
    - subsensor_DG_r.m (p<=r)
    - subsensor_DG_p.m (p>r)
  - sensor_QD.m
  
  - calculation_error.m (calculating percent reconstruction error)
    - subcalculation_error.m
  - calculation_det.m
  - averaged_operation.m
    - save_operation.m
  - data_organization.m
  
  - makevideo_NOAASST.m (making video)
  - makevideo_ensotrue.m (making video)
    - video_maker.m
      - display_sensors_color_for_movie.m
      
## How to cite
If you use the Determinant-based Fast Greedy Sensor Selection Algorithm code in your work, please cite the software itself and relevent paper.
### General software reference:
```bibtex
@misc{saito2019github,
      author = {Yuji Saito},
      title = {Determinant-based Fast Greedy Sensor Selection Algorithm},
      howpublished = {Available online},
      year = {2019},
      url = {https://github.com/YujiSaitoJapan/Determinant-based-Fast-Greedy-Sensor-Selection-Algorithm}
}
```
### Relevent paper reference:
```bibtex
@misc{saito2019determinantbased,
      title={Determinant-based Fast Greedy Sensor Selection Algorithm}, 
      author={Yuji Saito and Taku Nonomura and Keigo Yamada and Keisuke Asai and Yasuo Sasaki and Daisuke Tsubakino},
      year={2019},
      eprint={1911.08757},
      archivePrefix={arXiv},
      primaryClass={eess.SP}
}
```
## Author
SAITO Yuji

[Experimental Aerodynamics Laboratory](http://www.aero.mech.tohoku.ac.jp/eng/)
Department of Aerospace Engineering
Graduate School of Engineering

Tohoku University, Sendai, JAPAN

E-mail: yuji.saito@tohoku.ac.jp

Github: [YujiSaitoJapan](https://github.com/YujiSaitoJapan)
