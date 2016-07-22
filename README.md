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

- [1] [_The Kirby Institute_](https://kirby.unsw.edu.au/), The University of New South Wales, Sydney NSW 2052, Australia.
- [2] Epidemiology Unit VIDRL, The Doherty Institute, Melbourne Vic 3000, Australia. 
- [3] Department of Medicine, University of Melbourne, Victoria.

### Project organization ###

All the project files are stored in the main directory and 6 main sub-directories. This project has two layers of organization with the main model files contained in this parent directory as R markdown files. Files related to specific projects are contained within the projects directory.

_Main directory scripts_

The following R markdown scripts create specific projects, run the HBV model, and produce outputs. The numbering is to indicate the order the scripts need to be run. Each of the scripts require user inputs to specify project details. 

- **0-CreateProject.Rmd**: This script is used to setup a project based on the specifications entered by the user. The user has specify a project name, a start and end year for simulations, the simulation time step, and population names (which determines the number of populations in the model). The script creates a project directory in the "projects" directory and creates a project specifications .rda file and relevant sub-directories. Template input .csv files are also created for the user to input parameter values:
	- initial\_populations.csv: For specifying initial population sizes.
	- parameters\_constants.csv: For specifying parameters which constant over time.
	- parameters\_time\_varying.csv: For specifying the base value of parameters which change over time.
	- time\_varying\_ranges.csv: For specifying the relative uncertainty in the time varying parameters. This is done by specifying a range form each parameter at the start and end simulation times. The model does a linear interpolation between values sampled from these ranges. 
	- population\_interactions.csv: Specifies how populations interact with each other. If population X interacts with population Y then a 1 is put in the cell corresponding to row X and column Y. Generally this is symmetric but does not have to be. 
	- population\_transitions.csv: Specifies the rate populations moves from one to another (e.g. through aging). If population X moves into population Y at rate r then r is entered in the cell corresponding to row X and column Y.
	
- **1-SetupModel.Rmd**: This script is used to read in all the specified projects input .csv files and create the parameter sets that will be run in the model. This script needs to be run for the actual model simulations to be performed. The project parameter sets consist of a "best estimates" parameter set and an ensemble of parameter sets generated through sampling from the input parameter ranges. All parameter sets are appended to the project's .rda file. If the user changes or updates the .csv input files then this script need to be rerun for the correct parmaters to be used. 

- **2-RunHBVmodel.Rmd**: This script is used to actually run the HBV model for the specified project. All simulation results are appended to the project's .rda file.

- **3-SummaryResults.Rmd**: This script is used to produce figures of summary results for a specified project. It is written as an R markdown document so the figures can be generated as part of a summary results Word document using knitr. 

_Main directory sub-directories_

#### code ####

Contains specific functions and scripts used by the main modeling scripts. This directory also contains the original files for Ben Cowie's Berkeley Madonna model. 

#### data ####

Contains the data sets used in the original Berkeley Madonna model as well as any raw data relevant to the project. 

#### docs ####

Contains any documents we have written relation to the model and its development such as manuscripts, papers, reports, etc. Project specific documents are contained within the associated projects directory.

#### projects ####

The projects directory contains a sub-directory for each individual project created and run by the model. Each project directory in turn has the following structure:

- .csv files containing all the inputs to the model for the specific project.
- a project .rda file which contains all the project specifications, parameters, and results for the specific project.
- a results directory where project specific figures and outputs are stored.
- a docs directory where project specific documents (as per the main docs directory) are stored.
- a data directory where project specific data to inform input parameters or for calibration is stored.






