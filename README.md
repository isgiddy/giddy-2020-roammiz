# giddy-2020-roammiz
Analysis followed by Giddy et al 2020 - stirring by sea ice-melt water enhances submesoscale fronts in the Southern Ocean

At the moment, to implement the code here you need to clone this repository to your local machine.  
'git clone https://github.com/isgiddy/giddy-2020-roammiz.git'  

The data can be accessed via ftp. In your terminal,  
ftp ssh.roammiz.com  
Name: anonymous  

No password is required.

The relevant folder is giddy_2020
  
Load this data into the same folder as the repository, ensuring that the folder is named 'data'  

A python environment is provided which has all the relevant packages to enable you to run the code smoothly. Note this environment has only been tested in OS systems  

Not all the code is written by myself. I have relied heavily on the [GliderTools](https://pypi.org/project/glidertools/), spectra calculations are computed using code provided by J. DeShayes.  
   
Disclaimer: Although I have tried to make the code as generic as possible,
it is specific to the analysis followed by Giddy-etal-2020. For other implementations you will probably need to make a number of modifications.  
I hope the annotations within the code assist with this. 
