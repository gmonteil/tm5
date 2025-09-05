# Tasks (by order of priority)

- GUI (setup)
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
    - [ ] (1) GM by Wednesday to replace continuous line for obspack by dots or crosses in the GUI and other small things to improve visibitiy
- GUI (analysis)
    - [ ] (8) GM: check domain size on plot after change of region
    - [ ] (2) GM by Wednesday: Automated observation fit statistics (RMSE, maybe as table), first for pre-computed results, does the statistics adapt when we zoom on the time axis?
    - [ ] (3) GM after Wednesday: Enable selection between multiple sensitivity experiments (via a rearrangement of current GUI elements)
    - [ ] (7) GM: Integrate some visualisation of the emission budget in configuration part of the GUI (to allow perform a quick sanity check)
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations
    - [ ] improve the location and extend the content of the config file  
    - [ ] extend station explorer so the user can select time series from a list of completed experiments
    - [ ] implement submitting with same job ID would overwrite the previous job with the same ID, to avoid submitting the same job twice, to save resources (not needed immediately)
    - [ ] think about permission for co-working on same code base
- TM5 speedup
    - [ ] GM: speedup reading emissions
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] implement and test reading of the CAMS reanalysis OH field (to be used for the troposhere, for the stratosphere we'll rather use Bruehl et al.), this should be the default, Spivakovsky is kept-
    - [ ] (10) MV: check how Vaganovo looks now in global 1 by 1 run (launch on COSMOS), not too urgent
    - [ ] (1) GM to send to Arjo email with the macros with use, and get him to confirm 
- Sensitivity Runs (the important point is that we want to see an effect and the comparison with CAMS clearly shows us which run is more consistent)
    - [ ] (6) MV try cases with the regional emissions where in the sensitivity run we remove the emissions from:
                (c) a power plant near a site
- Documentation
    - [ ] properly setup github pages with mkdocs
- Inversion
    - [ ] implement dummy inversion based on Jacobians to test computational performance -> txk
    - [ ] (4) GM: adapt adjoint to run with the simplified emission preparation 
    - [ ] (5) GM: adapt adjoint to run with the OH field (to Spivakovsky first, which is enough for the Jacobian computations, possibly later CAMS OH) 
 
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- The combination of the GUI on the ICOS carbon portal with a different VM (like CODE.DE) should be simple to implement, it would require
  - migration of the VM part
  - mounting folder from additional VM on notebook
  - mounting folder from Jupyter hub on new VM not needed, because we would only need it to pass the yaml file to the VM, and that can be arranged differently ...

