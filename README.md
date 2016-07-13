## HBV Estimates and Projections Model ##

This project involves the development of an epidemiological model of hepatitis B virus (HBV) transmission in Australia. The primary purpose of the model is to provide estimates for the prevalence, incidence, and disease burden of HBV in Australia nationally, by each state and territory, and by local health district (LHD) for NSW and project future epidemic trajectories for each region. 
The project involves the collaboration between the Kirby Institute, UNSW Australia, and the Cowie group at the Epidemiology Unit VIDRL, The Doherty Institute, Victoria. NSW Health provided funding as part of the [_BBV & STI Research, Intervention and Strategic Evaluation Program (BRISE)_](http://www.brise.com.au/). While the model has been primarily developed to estimate the number of people living with chronic HBV and assocaited disease burden, the code has been designed to be as flexible as possible to faciliate future model development and applications. 

The HBV model is based on a previous model implemented in [_Berkeley Madonna_](http://www.berkeleymadonna.com/) and has been published previously (MacLachlan et. al. The Burden of Chronic Hepatitis B Virus Infection in Australia. Aust NZ J Public Health 37; 5: 416–22. doi:10.1111/1753-6405.12049). The current model is written in the R language as scripts or R markdown files and is currently **under development**.

**Corresponding author:** Richard T. Gray

**ORCID ID:** orcid.org/0000-0002-2885-0483

**Affiliation:** [_The Kirby Institute_](https://kirby.unsw.edu.au/), UNSW Australia, Sydney NSW 2052, Australia

### Contributors ###

The following people were the main contributors to the development of the HBV transmission model. 

- Dr Richard T. Gray: Principle investigator for the project at the Kirby Institute. Model developer and coder and maintainer of this repository [1].
- Assoc. Prof. Benjamin C. Cowie: Original developer of the model this code is based on. Principle investigator for the Cowie group at the Doherty Institute [2,3].
- Neil Bretana: Model developer and coder [1].
- Dr Jennifer H. MacLachlan: Model developer [2,3].
- Dr Nicole L. Allard: Model developer [2,3].

**Affiliations**:

-[1] [_The Kirby Institute_](https://kirby.unsw.edu.au/), The University of New South Wales, Sydney NSW 2052, Australia
-[2] Epidemiology Unit VIDRL, The Doherty Institute, Melbourne Vic 3000, Australia 
-[3] Department of Medicine, University of Melbourne, Victoria

### Project organization ###

All the project files are stored in 6 main directories. This project has two layers of organization with the main model files contained in this parent directory as R markdown files. Files related to specific projects are contained within the projects folder.

#### code ####

Contains specific functions and scripts used by the main model. This folder also contains the original files for Ben Cowie's Berkeley Madonna model. 

#### data ####

Contains the data sets used in the original Berkeley Madonna model as well as any raw data relevant to the project. 

#### docs ####

- manuscripts - papers, reports, etc based on the results of the model
- Other documents we are writing. 

#### misc ####

Contains any miscellaneous files associated 

#### outputs ####

#### projects ####



