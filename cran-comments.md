## R CMD check results

0 errors | 0 warnings | 3 notes

* checking installed package size ... NOTE
    installed size is 6.9Mb
    sub-directories of 1Mb or more:
      data   6.4Mb
      
  This is due to the included dataset necessary for package functionality.

* checking for future file timestamps ... NOTE
  unable to verify current time
  
  This appears to be a system-specific check issue.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'BayesMendel Lab <nkubista@hsph.harvard.edu>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Penetrance (3:8)
    penetrance (11:46)
    
  These are correctly spelled technical terms in genetics.

* The Title field starts with the package name:
  - This is intentional as the package name describes its primary function. 