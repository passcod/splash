# Splash

**[SPLAT!]** is a wonderful tool for RF analysis over terrain. Unfortunately,
it only works with a few terrain sources and only uses one core for processing.

[SPLAT!]: http://www.qsl.net/kd2bd/splat.html

**Splash** is three things, in this order:

1. A learning project to teach myself about the maths and physics underlying
   the [Longley-Rice] Irregular Terrain Model (ITM) as I go about porting it to
   Rust.
2. A learning project to teach myself deep optimisation over modern hardware,
   using at least SIMD and hopefully GPU to dramatically improve modelisation
   performance.
3. A production-grade interface to the ITM, and eventually other models.

[Longley-Rice]: https://en.wikipedia.org/wiki/Longley%E2%80%93Rice_model

## Porting goals

 - Not merely a port of the C++ port of the FORTRAN original, but a full
   refactor, with all functions and variables named sensibly.
 - Everything thoroughly documented inline, with the minor goal that one should
   be able to understand (at a high level) how the ITM works “simply” by
   reading Splash’s source.
 - 1950-era approximations replaced by exact versions for math functions where
   possible.
 - Entirely safe code.
 - Separation of preparation and execution routines.

## Research and papers

### Longley-Rice ITM

 - [The Irregular Terrain Model]: FORTRAN implementation with descriptions (NTIA 2002)
 - [Hufford 1995]: The (ITM) Algorithm, v1.2.2
 - [Hufford-Longley-Kissick 1982]: A Guide to the Use of the ITM
 - [McKenna 2016]: Propagation Measurement Workshop (Slides)

### Background

 - [Norton 1959]: Transmission Loss in Radio Propagation II
 - [Handbook on Ground Wave Propagation]: Edition of 2014 (ITU Radiocommunication Bureau)
 - [Phillips 2012]: Geostatistical Techniques for Practical Wireless Network Coverage Mapping

### Parallel ITM

 - [Song 2011]: Parallel Implementation of the Irregular Terrain Model (ITM) for Radio Transmission Loss Prediction Using GPU and Cell BE Processors
 - [Musselman 2013]: An OpenCL Implementation of the Irregular Terrain with Obstructions Model (ITWOM)

### Accuracy studies

 - [Sun 2015]: Propagation Path Loss Models for 5G Urban Micro-and Macro-Cellular Scenarios
 - [Sun 2016]: Investigation of Prediction Accuracy, Sensitivity, and Parameter Stability of Large-Scale Propagation Path Loss Models for 5G Wireless Communications
 - [Abhayawardhana 2005]: Comparison of Empirical Propagation Path Loss Models for Fixed Wireless Access Systems

### Other models or improvements

 - [El-Sallabi 2011]: Terrain Partial Obstruction LOS Path Loss Model for Rural Environments
 - [Phillips 2012]: Geostatistical Techniques for Practical Wireless Network Coverage Mapping
 - [MacCartney 2017]: Rural Macrocell Path Loss Models for Millimeter Wave Wireless Communications

[The Irregular Terrain Model]: https://www.its.bldrdoc.gov/media/50674/itm.pdf
[Hufford 1995]: https://www.its.bldrdoc.gov/media/50676/itm_alg.pdf
[Hufford-Longley-Kissick 1982]: https://www.ntia.doc.gov/files/ntia/publications/ntia_82-100_20121129145031_555510.pdf
[McKenna 2016]: https://www.its.bldrdoc.gov/resources/workshops/propagation-measurement-workshops-webinars.aspx
[Norton 1959]: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote12.pdf
[Handbook on Ground Wave Propagation]: https://www.itu.int/dms_pub/itu-r/opb/hdb/R-HDB-59-2014-PDF-E.pdf
[Phillips 2012]: https://core.ac.uk/display/54849067
[Song 2011]: https://ieeexplore.ieee.org/document/5680900/
[Musselman 2013]: https://github.com/amusselm/Parallel-LRP/blob/master/documents/report.pdf
[Sun 2016]: https://arxiv.org/abs/1603.04404
[Sun 2015]: https://arxiv.org/abs/1511.07311
[Abhayawardhana 2005]: https://ieeexplore.ieee.org/document/1543252
[El-Sallabi 2011]: https://ieeexplore.ieee.org/document/5701648
[MacCartney 2017]: https://ieeexplore.ieee.org/document/7914696
