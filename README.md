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
---
### Main program  
- P_greedy_demo.m  

### Function  
#### Preprocessing  
- F_pre_read_NOAA_SST.m  
- F_pre_SVD_NOAA_SST.m  
- F_pre_truncatedSVD.m  

#### Sensor selection  
- F_sensor_random.m  
- F_sensor_DC.m  
  - F_sensor_DC_sub.m  
    - F_sensor_DC_approxnt_vec.m  
    - F_sensor_DC_approxnt.m  
    - F_sensor_DC_loc.m  
    - F_sensor_DC_locr.m  
- F_sensor_QR.m  
  - F_sensor_QR_pivot.m  
- F_sensor_DG.m  
  - F_sensor_DG_r.m  
  - F_sensor_DG_p.m  
- F_sensor_QD.m  

#### Calculation
- F_calc_det.m  
- F_calc_sensormatrix.  
- F_calc_error.m  
  - F_calc_reconst.m  
  - F_calc_reconst_error.m  

#### Data organization  
- F_data_ave1.m  
- F_data_ave2.m  
- F_data_ave3.m  
- F_data_arrange1.m  
- F_data_arrange2.m  
- F_data_arrange3.m  
- F_data_arrange4.m
- F_data_arrange5.m
- F_data_normalize.m  

#### Mapping
- F_map_original.m  
	- F_map_videowriter.m  
		- F_map_plot_sensors_forvideo.m  
- F_map_reconst.m  
	- F_map_plot_sensors.m  

### Function  
#### Preprocessing  
- F_pre_read_NOAA_SST.m  
- F_pre_SVD_NOAA_SST.m  
- F_pre_truncatedSVD.m  

      
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
