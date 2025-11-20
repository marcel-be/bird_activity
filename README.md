## Workflow to analyse (movement) activity pattern of VHF tacked animals.

Here I provide R Scripts for a study on a forest bird community. For more information about the underlying study and data curation read *Lindner & Becker et al. (in prep.)*.
In our study, birds were equipped with VHF transmitters, the data wereThe data were recorded and processed via the tRackIT system (see *Höchst et al. 2021* for further details). 

The basis for this script is a dataset with one-minute activity bouts for each bird individual and date, with each minute labelled as either ‘active’ or ‘passive.’ (see *Gottwald et al. 2022* for further details).
For each bird, daily activity probabilities were estimated by fitting hierarchical GAMs, generating individual-level and day-specific ‘activity curves.’. From the predicted activity curves, we then derived certain activity characteristics for each day of each individual. 

To quantify how variation in daily activity characteristics is partitioned across biological levels-between species, between individuals within species, and within individuals across days-we conducted a variance-components analysis based on hierarchical linear mixed-effects modelling

In addition to all analyses, I provide code for the graphical evaluation of the activity curves and diagrams

### Literature for further details

Gottwald, J., Royauté, R., Becker, M.,..., 2022. Classifying the activity states of small vertebrates using automated VHF telemetry. https://doi.org/10.1101/2022.03.22.485147

Höchst, J., Gottwald, J., Lampe, P., ..., 2021. tRackIT OS: Open-source Software for Reliable VHF Wildlife Tracking. Presented at the INFORMATIK 2021, Gesellschaft für Informatik, Bonn, pp. 425–442.
