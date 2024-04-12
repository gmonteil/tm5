I put here some of the files from the base that get removed from the source by pycasso.

- Files directly in this directory were removed by "build.configure.remove"
- Files in subdirectories were removed by macros (build.configure.macro.*), e.g. *dry\_deposition* contains the files that used to be removed by `without_dry_deposition`.
- Files in the _skip_ subdirectory are those that get removed by the "build.copy.skip.file" rc-key

