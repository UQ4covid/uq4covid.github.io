The script files and data sources used here are to inform age-specific contact matrices.

Data (`20200327_comix_social_contacts.xlsx`) from here:

(https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html)[https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html]

Setting the working directory to this folder, and running the
`contactData.R` script produces a file called `coMix_matrix.csv`
for use in the model which corresponds to the 
`All_contacts_imputed` sheet with column `A` and row `1` removed.

We also extract the POLYMOD survey data form the `socialmixr` package,
and save this as the file `POLYMOD_matrix.csv`.
