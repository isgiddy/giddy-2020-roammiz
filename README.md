# giddy-2020-roammiz
Analysis followed by Giddy et al 2020 - stirring by sea ice-melt water enhances submesoscale fronts in the Southern Ocean

[![DOI](https://zenodo.org/badge/288300583.svg)](https://zenodo.org/badge/latestdoi/288300583)

At the moment, to implement the code here you need to clone this repository to your local machine.  
'git clone https://github.com/isgiddy/giddy-2020-roammiz.git'  

The data can be accessed via ftp. In your terminal,  
ftp ssh.roammiz.com  
Name: anonymous  

No password is required.

The relevant folder is giddy_2020
  
Load this data into the same folder as the repository, ensuring that the folder is named 'data'  

A python environment is provided which has all the relevant packages to enable you to run the code smoothly. Note this environment has only been tested in OS systems  

Thanks go out specifically to the Xarray community for making life so much easier (and apologies for not using Xarray optimally) and generally to the open-code community.  
For reading Seaglider data, I have relied heavily on [GliderTools](https://pypi.org/project/glidertools/), spectra calculations are computed using code provided by J. DeShayes.  
   
Disclaimer: The purpose of this repository is to make public the analysis followed in Giddy-etal-2020 under the philosophy of open-code, reproducibility and transparency. Although I have tried to make the code as generic as possible,
it is specific to the analysis followed by Giddy-etal-2020. For other implementations you will probably need to make a number of modifications.  
I hope the annotations within the code assist with this. 
