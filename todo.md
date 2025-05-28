# Tasks (by order of priority)

- GUI (setup)
    - [x] Test on ICOS Jupyter hub and fix problem
    - [x] Make conda environment accessible for the FIT-IC group
    - [x] Ask inventory compilers to register at Jupyter hub and ask Ute to make them members of fit-ic project
    - [x] (re) implement connection with TM5 (work with Zois)
    - [x] add GUI code to TM5 code on this GIT
    - [x] add an "overwrite" field for the emissions
    - [ ] change: "Use different emissions for the zoom domain" -> "Use different regional emissions"
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
- GUI (analysis)
    - [x] plot continuous timeseries at requested stations
    - [ ] check domain size on plot after change of region
    - [ ] Automated observation fit statistics (RMSE, maybe as table)
    - [x] automated plot with emissions budget by region and category
    - [ ] Integrate some visualisation of the emission budget in configuration part of the GUI (to allow perform a quick sanity check)
    - [x] plot comparison with observations
- TM5
    - [x] diagnose issue with concentrations decreasing too fast (i.e. sink too strong?)
    - [x] -> mv to check preprocessing of emission for climatological files
    - [x] -> mv to test 1 year simulation with regional zoom (takes roughly 1 hour on pancake and on COSMOS, with single tracer)
    - [x] -> mv to prepare the met forcing for the zoom configurations
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations
    - [x] implement a job scheduler so multiple TM runs triggered by different users can happen at the same time
    - [ ] improve the location and extend the content of the config file  
    - [ ] extend station explorer so the user can select time series from a list of completed experiments
    - [ ] implement submitting with same job ID would overwrite the previous job with the same ID, to avoid submitting the same job twice, to save resources (not needed immediately)
    - [ ] think about permission for co-working on same code base
- TM5 speedup
    - [x] speedup chemistry
    - [ ] speedup reading emissions
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] check proper reading of the OH file and see whether that improves the seasonality at SPO -> GM
    - [ ] implement and test reading of the CAMS reanalysis OH field (to be used for the troposhere, for the stratosphere we'll rather use Bruehl et al.), this should be the default, Spivakovsky is kept as an alternative -> GM
    - [ ] simulated variability of Cabauw appears too small; do forward run with out model setup and posterior emission fields from CAMS inversion -> MV to check for the location of the emissions; GM to check in his run with Avengers zoom configuration
- Documentation
    - [ ] properly setup github pages with mkdocs
 
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- The combination of the GUI on the ICOS carbon portal with a different VM (like CODE.DE) should be simple to implement, it would require
  - migration of the VM part
  - mounting folder from additional VM on notebook
  - mounting folder from Jupyter hub on new VM not needed, because we would only need it to pass the yaml file to the VM, and that can be arranged differently ...

