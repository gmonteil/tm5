# Tasks (by order of priority)

- prior
    - [x] add comment or extend long name forfield "fossil" in EDGAR files to clarify that this is the sum of a number of fields in the same file(s)
- GUI (setup)
    - (6) [x] MV: change: "Use different emissions for the zoom domain" -> "Use different regional emissions"
    - [ ] implement only configuration with 3 level zoom "avengers-1", i.e global/Europe/NL+D+CH
    - [2] GM to replace continuous line for obspack by dots or crosses 
    - [3] GM to check again unit conversion for zoom region emissions
    - [4] MV to check the preparation of the monthly emissions files (monthly and flat) 
- GUI (analysis)
    - (8) [ ] GM: check domain size on plot after change of region
    - [ ] Automated observation fit statistics (RMSE, maybe as table)
    - [x] automated plot with emissions budget by region and category
    - (7) [ ] GM: Integrate some visualisation of the emission budget in configuration part of the GUI (to allow perform a quick sanity check)
    - [x] plot comparison with observations
- GUI (extension)    
    - [ ] (re) implement possibility to pass defaults through a YAML file
    - [ ] implement a few default configurations
    - [ ] improve the location and extend the content of the config file  
    - [ ] extend station explorer so the user can select time series from a list of completed experiments
    - [ ] implement submitting with same job ID would overwrite the previous job with the same ID, to avoid submitting the same job twice, to save resources (not needed immediately)
    - [ ] think about permission for co-working on same code base
- TM5 speedup
    - (5) [ ] GM: speedup reading emissions, and check correct implementation in zoom; also check code against original version from https://sourceforge.net/p/tm5/cy3_4dvar/ci/cams-ch4/tree/ (if this does not work then we interact more with Arjo)
- TM5 analyses, open questions regarding our initial TM5 runs with prior emission fields: to be addressed with long run on COSMOS:
    - [ ] implement and test reading of the CAMS reanalysis OH field (to be used for the troposhere, for the stratosphere we'll rather use Bruehl et al.), this should be the default, Spivakovsky is kept-
    - (1) [ ] MV: make comparison between 3 pm obs at Cabauw and simulated value at same time of day
- Documentation
    - [ ] properly setup github pages with mkdocs
- Inversion
    - [ ] implement dummy inversion based on Jacobians to test computational performance -> txk
    - [ ] adapt adjoint to run with the simplified emission preparation -> gm
    - [ ] adapt adjoint to run with the OH field (to Spivakovsky first, which is enough for the Jacobian computations, possibly later CAMS OH) -> gm
 
# Task ideas / discussion
- avoid misuse through implementation of "accepted ranges" => *That's essentially built-in the GUI*
- Documentation: could be overleaf or a shared word document ... https://lumia.nateko.lu.se/ using a package mkdoc, could be combined with a markedown file on nextcloud ... => *actually, there's already a squeletton in place. Just install the `mkdocs` python package and run `mkdocs serve` from within the root of the repo. The content of the resulting page is based on the markdown files in the `docs` folder*
- The combination of the GUI on the ICOS carbon portal with a different VM (like CODE.DE) should be simple to implement, it would require
  - migration of the VM part
  - mounting folder from additional VM on notebook
  - mounting folder from Jupyter hub on new VM not needed, because we would only need it to pass the yaml file to the VM, and that can be arranged differently ...

