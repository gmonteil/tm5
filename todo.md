# Tasks (by order of priority)

- GUI (setup)
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
- GUI (analysis)
    - [ ] (8) GM: check domain size on plot after change of region
    - [ ] (4) GM+MV: extract TM5 output and obspack for selected stations averaged over a pre-specified period of the day (as in inversions), e.g. above 1000 m height take nighttime data (midnight to 4am local solar time) and below afternoon date (from 12 to 4 local solar time) (there is a flag in obspack for filtering, but don't take it as an absolute reference, just for a sanity check) at highest level
    - [x partly] now: Automated observation fit statistics (RMSE, maybe as table with both RMSEs and a map also indicating the relative difference of two RMSEs), first for pre-computed results, does the statistics adapt when we zoom on the time axis? Add Chi2 with sigma = sqrt (sigma_obs**2 + sigma_model **2), and sigma model from a crude approximation derived from short scale variability of the obs, or the variability of the meteo (representation error approximation from slopes may also be there ...
    - [ ] (7) GM: Integrate some visualisation of the emission budget in configuration part of the GUI (to allow perform a quick sanity check)
    - [ ] GM: Add plot or table that indicates relative contribution of each sector to total signal at each site (only if we run sectors separately), this is not urgent, can be added later ...
    - [X] (3) GM+MV: Inspect what makes the first cell in the GUI so slow, must be some data processing that we should be beforehand, because it is too slow ...
    - [ ] (1) MV: More text around plots to explain what we see, e.g. time series show hourly values but fit statistics are based on daily averages (as a preparation for the selection of averges over specific periods of the day)
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations, add the option to transport emission categories separately (with a warning on performance)
    - [ ] improve the location and extend the content of the config file  
    - [ ] extend station explorer so the user can select time series from a list of completed experiments
    - [ ] implement submitting with same job ID would overwrite the previous job with the same ID, to avoid submitting the same job twice, to save resources (not needed immediately)
    - [ ] think about permission for co-working on same code base
- TM5 speedup
    - [ ] GM: speedup reading emissions
    - [ ] GM: 7200 seconds for 7 tracer run was much longer than for the single-tracer run (~4800 seconds)
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] implement and test reading of the CAMS reanalysis OH field (to be used for the troposhere, for the stratosphere we'll rather use Bruehl et al.), this should be the default, Spivakovsky is kept-
    - [ ] (10) MV: check how Vaganovo looks now in global 1 by 1 run (launch on COSMOS), not too urgent
- Sensitivity Runs (the important point is that we want to see an effect and the comparison with CAMS clearly shows us which run is more consistent)
    - [ ] (6) MV to further analyse the default vs regional cases by comparing signal from component fluxes offline (and without observations) and show these plots together with the plot using the sectoral totals and the observations that we already have in the tool
- Documentation
    - [ ] properly setup github pages with mkdocs
- Inversion
    - [ ] TXK: implement dummy inversion based on Jacobians to test computational performance
    - [ ] TXK: to specify how the footprints wil be used in the inverse model (dc/df, no emission caterogies).
    - [ ] TXK/MV: draft a kind of tutorial guiding the user through one of the experiments
    - [ ] (2a) GM: check units of footprint 
    - [ ] (2b) MV to further test footprint calculation
    - [ ] GM: make sure that adjoint file is returned in the right place
          
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- The combination of the GUI on the ICOS carbon portal with a different VM (like CODE.DE) should be simple to implement, it would require
  - migration of the VM part
  - mounting folder from additional VM on notebook
  - mounting folder from Jupyter hub on new VM not needed, because we would only need it to pass the yaml file to the VM, and that can be arranged differently ...
- Note that the concentration calculation uses fine grid in the full zoom domain, but the emissions close to the boundary are using the coarse resolution (that is applied on the other side of the boundary). If we wanted to change that we would need to enlarge the zoom domain.

