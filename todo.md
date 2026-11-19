# Tasks (by order of priority)
- GUI development steps
    - [ ] (1) There is a problem with the file upload from Notebooks on ICOS hub: When the file size passes a certain threshold (roughly 15MB), then it is not uploaded and no error message displayed. This is a Jupyter feature! A potential fix for us would be to use the upload tool separately and then access the uploaded file from our Notebook. But for now we can live with the restriction.
    - [ ] (2) MV/txk to adapt fortran code such that both the short demo version and the full one year version both rely on the same file for the Jacobian and the observations, the demo version just reads in and uses only the fraction needed to cover Jan 2021 for the obs and the period from Oct 2020 to Jan 2021 for the Jacobian (need to observe CPU time for the demo version, it should remain as quick as it is now)
    - [ ] (3) MV to make this line invisible to the user: "yamlfile='/srv/data/AVENGERS/fit_ic/gui/gui-new_ilabdesk.yml'"
    - [ ] (4) MV to adapt the GUI such that the choice of the configuration is between "Fast version Jan 2021" and "Slower version for full year 2021"
    - [ ] (5) MV/txk to adapt fortran code such that the choice under item (4) is read in and used
- GUI (analysis)
    - []  TK to come up with a strategy for this: Add Chi2 with sigma = sqrt (sigma_obs**2 + sigma_model **2), and sigma model from a crude approximation derived from short scale variability of the obs, or the variability of the meteo (representation error approximation from slopes may also be there ...
    - [ ] (1) GM: big cleanup!
    - [ ] (2) GM: Integrate some visualisation of the emission budget in configuration part of the GUI (to allow perform a quick sanity check)
    - [ ] (9) GM: to implement that option of making an uploaded emission file public (e.g. only show the user specific subdir in the emissions folder to that very user); maybe we also need to validate the file ourselves before it is allowed to be public
    - [ ] (10) GM: Add plot or table that indicates relative contribution of each sector to total signal at each site (only if we run sectors separately), this is not urgent, can be added later ...
    - [ ] The 'coarsen_file' routine for the emissions expected dimensions 'lon'/'lat' in the file, while the upload form introduced 'longitude'/'latitude'.
For now I have adapted this in the  upload_form.py, but eventually we probably should relax the naming conventions a bit. (not urgent)
- GUI (extension)
    - [ ] try to extend the current logic that allows us to retrieve results from previous experiments beyond the current session (use checksum of config file)
TM5:
    - [ ] (10) MV: check how Vaganovo looks now in global 1 by 1 run (launch on COSMOS), not too urgent
- Documentation
    - [ ] properly setup github pages with mkdocs
- Inversion
    - [ ] MV to extract the wrong and the corrected version of the TM5 subroutine (reconstruct_required_var_1) that was responsible for the error in the adjoint
    - [ ] MV to test GM's update of stratospheric chemistry     
    - [ ] TXK/MV: draft a kind of tutorial guiding the user through one of the experiments
    - [ ] TXK to exclude spinup from target Jacobian
    - [ ] When we run the footprints for obs from Jan to Dec 2021 back to October 2020 emissions we do the following: For filtering of the footprints wrt to the network, we include the information about the the network into the observation table file  
            
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- The combination of the GUI on the ICOS carbon portal with a different VM (like CODE.DE) should be simple to implement, it would require
  - migration of the VM part
  - mounting folder from additional VM on notebook
  - mounting folder from Jupyter hub on new VM not needed, because we would only need it to pass the yaml file to the VM, and that can be arranged differently ...
- Note that the concentration calculation uses fine grid in the full zoom domain, but the emissions close to the boundary are using the coarse resolution (that is applied on the other side of the boundary). If we wanted to change that we would need to enlarge the zoom domain.
- For the global Jacobian we use the NOAA flask sampling network, but we compute the footprints for all flasks in Obspack.

