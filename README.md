# RELAGN
Spectral model for the calculation of AGN SEDs, ranging from the Optical/UV (outer accretion disc) to the Hard X-ray (Innermost X-ray Corona), including
a fully relativistic treatment. If this code is useful to your work, then please reference: **Hagen &amp; Done (2023b, submitted)**

This model is based off the physically motivated broad-band SED model AGNSED (Kubota &amp; Done 2018), however with the addition of a fully relativistic
treatment. The emissivity of the accretion flow is calculated following the relativistic treatment of Novikov &amp; Thorne (1973) (This is also done in AGNSED),
while the ray tracing from the point of emission to the observer is performed using the relativistic convolution code KYCONV (Dociak, Karas &amp; Yaqoob 2004.)
For a detailed description of the model please see Hagen & Done (2023b, submitted) (https://ui.adsabs.harvard.edu/abs/2023arXiv230401253H/abstract)

There are currently two version of the code: one in Python and another in Fortran. The Fortran version is written to be used with the spectral fitting
software XSPEC. If you intend to use the model to analyse X-ray spectral data, then this will probably be the prefered version to use. The Python version,
on the other hand, provides more flexibility for modelling purposes. Where the Fortran version will only produce a spectrum, the python version can be used
to easily extract the physical properties of the system (e.g physical mass accretion rate, disc size, efficiency parameters, etc..) since these are all stored
as attributes within the model (it's written in an object oriented manner). This would be more appropriate for working outside of XSPEC. 
Both versions have the same input parameters, and are detailed below.

It is also important to note, that both versions require a working installation of HEASOFT, since we use the XSPEC convolution model KYCONV to apply
the relativistic effects. If you do not already have HEASOFT installed, it can be found here: https://heasarc.gsfc.nasa.gov/docs/software/heasoft/

Requirements
------------
* Working installation of HEASOFT, since we use the XSPEC convolution model KYCONV to apply the relativistic effects this is required for both the
Python and Fortran versions. If you do not already have HEASOFT installed, it can be found here: https://heasarc.gsfc.nasa.gov/docs/software/heasoft/.
The code has been tested for HEASOFT versions: 6.29, 6.30, and 6.31; with XSPEC versions 12.12.0, 12.12.1, and 12.13.0 </br>
(If you have succesfully installed HEASOFT then you will have everything you need to run the Fortran version - seeing as HEASOFT won't build 
without the necessary compilers...)

* Numpy (ONLY for Python version - Tested for v.1.21.5)
* Scipy (ONLY for Python version - Tested for v.1.9.3)
* Astropy (ONLY for Python version - Tested for v.5.1)

Note: The Python version will ONLY run in Python 3! - We have been using v.3.9.12


Installation
-------------
### Fortran
1. Download the source code files `relagn.f` and `lmod_relagn.dat`, found within `~/src/fortran_version`
2. Place in a directory of your choice, and open XSPEC
3. Whthin XSPEC type: `initpackage relagn lmod_relagn.dat .` This will compile the code. Note, this assumes you are currently within the directory
containing the source code files...
4. Still within XSPEC type: `lmod relagn .` This will load RELAGN as a local model in XSPEC (call as relagn). Again, this assumes you are currently
within the directory containing the source code. If you wish to load it from a different location, replace `.` with the explicit path to the relagn 
directory (i.e `/path/to/relagn`)
5. (OPTIONAL) By default you will have to repeat step 4. every-time you open XSPEC. This is tedeous. Instead you can put the line `lmod relagn /path/to/relagn`
in your `xspec.rc` file; located within your `~/.xspec` directory. If `xspec.rc` does not exist within `~/.xspec`, you can simply make one!
6. (ALTERNATIVE) Instead of manully performing steps 3 to 5, you can simply excecute the shell scripts `compile_to_xspec.sh` and `init_autoLoad_xspec.sh`, which will perform steps 3-4 and 5 respectively


### Python
Currently the easiest way is to either clone this repository, or directly download the source code (i.e ALL the files within the `~/src/python_version`
directory). Then add the source code file to your PYTHONPATH (or use sys.append() within your own scripts if you do not wish to directly edit the
PYTHONPATH...). Turning this into a package you can pip install is on the to-do list!!!

Model Parameters
----------------
**Par 1. &ensp;  $M$** </br>
  &emsp; &emsp; &#9656; **Units:** $M\_{\odot}$ </br>
  &emsp; &emsp; &#9656; **Description:** Mass of the central Black Hole

**Par 2. &ensp;  $D$** </br>
  &emsp; &emsp; &#9656; **Units:** Mpc </br>
  &emsp; &emsp; &#9656; **Description:** Co-Moving distance from the observer to the Black Hole </br>
  
**Par 3. &ensp;  $\log \dot{m}$** </br>
  &emsp; &emsp; &#9656; **Units:** $\log \dot{M}/\dot{M}\_{\mathrm{Edd}}$ </br>
  &emsp; &emsp; &#9656; **Descripton:** Log Mass-accretion rate, scaled by the Eddington mass accretion-rate. 
  (i.e $\log \dot{m} = -1$ would imply $\dot{M} = 0.1 \dot{M}\_{\mathrm{Edd}}$)
 
 **Par 4. &ensp; $a$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description:** Black hole spin parameter. 0 implies non-spinning, while 1 is maximally spinning (with prograde rotation). 
  Note, the code will limit you to max 0.998 - This is the theoretical maximum assuming the presence of a disc

 **Par 5. &ensp; $\cos(i)$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description**: Cosine of the inclination of the observer with respect to the disc, as measured from the z-axis with 
  the disc in the x-y plane 
 
 **Par 6. &ensp; $kT\_{e, h}$** </br>
  &emsp; &emsp; &#9656; **Units:** keV </br>
  &emsp; &emsp; &#9656; **Description:** Electron temperature for the hot Comptonising corona. 
  This sets the high-energy roll-over for the hot Comptonisation region
 
 **Par 7. &ensp; $kT\_{e, w}$** </br>
  &emsp; &emsp; &#9656; **Units** keV </br>
  &emsp; &emsp; &#9656; **Description:** Electron temperature for the warm Comptonising region.
 
 **Par 8. &ensp; $\Gamma\_{h}$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description:** Spectral index for the hot Comptonisation component
 
 **Par 9. &ensp; $\Gamma\_{w}$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description:** Spectral index for the warm Comptonisation component
 
 **Par 10. &ensp; $r_{h}$** </br>
  &emsp; &emsp; &#9656; **Units:** $R\_{G}$ ( $R\_{G} = GM/c\^{2}$ so technically dimensionless) </br>
  &emsp; &emsp; &#9656; **Description:** Outer radius of the hot Corona. 
  If this is negative, then the code will set it to the innermost stable cirbular orbit, $r\_{\mathrm{isco}}$
 
 **Par 11. &ensp; $r_{w}$** </br>
  &emsp; &emsp; &#9656; **Units:** $R\_{G}$ </br>
  &emsp; &emsp; &#9656; **Description:** Outer radius of the warm Comptonisation region. 
  If this is negative, then the code will set it to $r\_{\mathrm{isco}}$
 
 **Par 12. &ensp; $\log r\_{\mathrm{out}}$** </br>
  &emsp; &emsp; &#9656; **Units:** $R\_{G}$ </br>
  &emsp; &emsp; &#9656; **Description:** Log of the outermost disc radius. 
  If this is negative, then the code will use the self-gravity radius from Laor &amp; Netzer (1989)
 
 **Par 13. &ensp; $f\_{\mathrm{col}}$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description:** Colour-temperature correction to be applied to the standard **outer* disc. 
  If this is negative, then the code follows the relation in Done et al. (2012). Otherwise treated as a constant correction
 
 **Par 14. &ensp; $h\_{\mathrm{max}}$** </br>
  &emsp; &emsp; &#9656; **Units:** $R\_{G}$ </br>
  &emsp; &emsp; &#9656; **Description:** Maximum scale-height of the inner-corona
 
 **Par 15. &ensp; $z$** </br>
  &emsp; &emsp; &#9656; **Units:** Dimensionless </br>
  &emsp; &emsp; &#9656; **Description:** Redshift of the source
 
 </br>
 
 For more details, see Hagen & Done (2023b, submitted)
 
Citing RELAGN
-------------
If you use RELAGN in your work, please cite Hagen & Done (2023b, submitted). The paper is currently in peer review, so for now please cite the version on the arXive (https://ui.adsabs.harvard.edu/abs/2023arXiv230401253H/abstract). You can use the following bibtex:

```
@ARTICLE{2023arXiv230401253H,
       author = {{Hagen}, Scott and {Done}, Chris},
        title = "{Estimating Black Hole Spin from AGN SED Fitting: The Impact of General-Relativistic Ray Tracing}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Astrophysics of Galaxies},
         year = 2023,
        month = apr,
          eid = {arXiv:2304.01253},
        pages = {arXiv:2304.01253},
          doi = {10.48550/arXiv.2304.01253},
archivePrefix = {arXiv},
       eprint = {2304.01253},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230401253H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

