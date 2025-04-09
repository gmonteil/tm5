# Tasks (by order of priority)

- GUI (setup)
    - [x] Test on ICOS Jupyter hub and fix problem
    - [x] Make conda environment accessible for the FIT-IC group
    - [x] Ask inventory compilers to register at Jupyter hub and ask Ute to make them members of fit-ic project
    - [ ] (re) implement connection with TM5 (work with Zois)
    - [x] add GUI code to TM5 code on this GIT
    - [ ] add an "overwrite" field for the emissions
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
- GUI (analysis)
    - [x] plot continuous timeseries at requested stations
    - [ ] check domain size on plot after change of region
    - [ ] Automated observation fit statistics (RMSE, maybe as table)
    - [ ] automated plot with emissions budget by region and category
    - [x] plot comparison with observations
- TM5
    - [x] diagnose issue with concentrations decreasing too fast (i.e. sink too strong?)
    - [x] -> mv to check preprocessing of emission for climatological files
    - [ ] -> mv to test 1 year simulation with regional zoom
    - [x] -> mv to prepare the met forcing for the zoom configurations
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations
    - [ ] implement a job scheduler so multiple TM runs triggered by different users can happen at the same time
- TM5 speedup
    - [x] speedup chemistry
    - [ ] speedup reading emissions
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] check proper reading of the OH file and see whether that improves the seasonality at SPO -> GM
    - [ ] simulated variability of Cabauw appears too small; do forward run with out model setup and posterior emission fields from CAMS inversion -> MV to check for the location of the emissions
- Documentation
    - [ ] properly setup github pages with mkdocs
 
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- Advanced GUI: Comparison between two (or more?) simulations
