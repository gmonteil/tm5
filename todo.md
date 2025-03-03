# Tasks (by order of priority)

- GUI (setup)
    - [x] Test on ICOS Jupyter hub and fix problem
    - [x] Make conda environment accessible for the FIT-IC group
    - [x] Ask inventory compilers to register at Jupyter hub and ask Ute to make them members of fit-ic project
    - [ ] (re) implement connection with TM5
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
- GUI (analysis)
    - [x] plot continuous timeseries at requested stations
    - [ ] plot comparison with observations
    - [ ] automated plot with emissions budget by region and category
    - [ ] Automated observation fit statistics
- TM5
    - [x] diagnose issue with concentrations decreasing too fast (i.e. sink too strong?)
    - [x] -> mv to check preprocessing of emission for climatological files
    - [ ] -> gm to test 1 month simulation with regional zoom
    - [x] -> mv to prepare the met forcing for the zoom configurations
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations
    - [ ] implement section(s) for observations
- TM5 speedup
    - [x] speedup chemistry
    - [ ] speedup emissions
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] repeat experiment with initial conditon, e.g. from CAMS, so the initial spatial gradient is realistic -> MV, txk to check 
    - [ ] simulated variability of Cabauw appears too small. It is primarily driven by changes in boundary layer height and wind direction -> MV, txk to check effect of vertical resolution
- Documentation
    - [ ] properly setup github pages with mkdocs
 
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- Advanced GUI: Comparison between two (or more?) simulations
