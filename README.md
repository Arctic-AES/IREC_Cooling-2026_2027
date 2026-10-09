# Steady-State Avionics Thermal 1D Model

Avionics-bay 1D thermal model for UGA AIAA's IREC rocket, Church.


| Board | Named components | Small parts | Total heat | Full-board IHG for supplied STEP body |
|---|---:|---:|---:|---:|
| C6 | 4.533 W | 0.194 W | **4.727 W** | **216,270 W/m³** |
| G6 | 0.600 W | 0.091 W | **0.691 W** | **31,587 W/m³** |
| ENIAC | 6.510 W | 0.133 W | **6.643 W** | **707,176 W/m³** |

Total: **12.061 W**. 

| Scenario | Air / ground | Exterior h | Solar | Bay air | C6 | G6 | ENIAC |
|---|---|---:|---|---:|---:|---:|---:|
| Hot pad | 42 / 60°C | 10 W/(m² K) | on | 67.56°C | 78.65°C | 66.97°C | 107.87°C |
| Moderate pad | 34 / 52°C | 21.6 W/(m² K) | on | 51.70°C | 64.06°C | 51.19°C | 95.84°C |
| Cool / no sun | 20.5 / 20.5°C | 22 W/(m² K) | off | 27.99°C | 42.55°C | 27.66°C | 78.52°C |


## Original reference list

The following references were carried forward from the original project README;
their inclusion does not certify every current numerical assumption.


1. ESRA, 2026 IREC Design, Test & Evaluation Guide V1.1, Appendix G. https://www.soundingrocket.org/what-is-irec.html
2. NOAA NCEI, Engineering Weather Data, Midland TX. https://www.ncei.noaa.gov/pub/data/EngineeringWeatherData_CDROM/engwx/midland_tx.pdf
3. Bergman, Lavine, Incropera, DeWitt, Fundamentals of Heat and Mass Transfer, 7th ed., Wiley: ch. 1, 3, 7 (Churchill-Bernstein), 9 (Churchill-Chu), 13, Table A.4.
4. PVEducation, Standard Solar Spectra (ASTM G173). https://www.pveducation.org/node/430
5. Wikipedia, Albedo. https://en.wikipedia.org/wiki/Albedo
6. Swinbank, W.C. (1963), Long-wave radiation from clear skies, Q. J. R. Meteorol. Soc. 89, 339-348. https://doi.org/10.1002/qj.49708938105
7. NASA POWER climatology, launch site 31.05 N, 103.55 W. https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=TS_MAX,T2M_MAX,ALLSKY_SRF_ALB&community=RE&longitude=-103.5473&latitude=31.0498&format=JSON
8. FLIR, What emissivity settings should I use when looking at PCBs and electronics? https://www.flir.com/en-ca/support/instruments2/what-emissivity-settings-should-i-use-when-looking-at-pcbs-and-electronics/
9. Wikipedia, FR-4. https://en.wikipedia.org/wiki/FR-4
10. Datasheets: TL780 https://www.ti.com/lit/ds/symlink/tl780.pdf, TPS563200 https://www.ti.com/lit/ds/symlink/tps563200.pdf, TLV1117LV https://www.ti.com/lit/pdf/SBVS160, LMX2541 https://www.ti.com/lit/ds/symlink/lmx2541.pdf, HMC361 https://www.analog.com/media/en/technical-documentation/data-sheets/hmc361s8g.pdf, BMT321 https://www.berex.com/Products/HighPowerAmplifier, SX1281 https://www.mouser.com/datasheet/2/761/DS_SX1280-1_V2.2-1511144.pdf, ADV7280A https://www.analog.com/media/en/technical-documentation/data-sheets/adv7280a.pdf, SiT5000 https://www.sitime.com/datasheet/SiT5000, NEO-M9N https://content.u-blox.com/sites/default/files/NEO-M9N-00B_DataSheet_UBX-19014285.pdf, MMC5983MA https://www.memsic.com/Public/Uploads/uploadfile/files/20220119/MMC5983MADatasheetRevA.pdf, ISL83488 https://www.renesas.com/us/en/document/dst/isl83483-isl83485-isl83488-isl83490-isl83491-datasheet
