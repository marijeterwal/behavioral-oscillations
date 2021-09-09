# behavioral-oscillations

This repository contains the <b>code</b> used for the manuscript:
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, R., Staresina, B., Hanslmayr, S., Wimber, M. 
Theta rhythmicity governs the timing of behavioural and hippocampal responses in humans specifically during memory- dependent tasks. 
bioRxiv 2020.11.09.374264; DOI: https://doi.org/10.1101/2020.11.09.374264

The <b>datasets</b> are available here:
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., 
Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, 
R., Staresina, B., Hanslmayr, S., Wimber, M., (2020), Data for: Theta 
rhythmicity governs the timing of behavioural and hippocampal responses in 
humans specifically during memory- dependent tasks. figshare. Collection. 
DOI: https://doi.org/10.6084/m9.figshare.c.5192567

The results were also presented as a poster at CNS 2020 (poster G183), titled 'Oscillatory patterns in behavioral and hippocampal responses during an associative memory task'. 
For a walk-through, see https://youtu.be/28kpbGDHLuo


## Prerequisites
MATLAB (The MathWorks) - The code has been tested on Matlab versions between R2018a and R2018b.

## Getting started
Download the .zip and unzip or clone to your favourite path.
Make sure you Matlab path is set to include the code.

Download the relevant data (see above) and adjust the paths in the scripts to point to the data.
You are ready to run the code! <br>
For more detailed instructions please see the analysis scripts in the folder Scripts.

The behavioral analyses should, on a regular computer, only take a few minutes to run per task phase. The PPC analyses and simulations (with the current settings) take more data/computational steps and will therefore take substantially more time (~hour(s)), depending on the available hardware.

You can of course also use the O-score function on your own data. Any time series can be entered into the O-score function. For details on how to set up the configuration please refer to the file oscillationScore.m.

## License and disclaimer
This code is released under a CC BY 4.0 license, meaning that this code is free to use and modify for everyone, given appropriate credit, but comes without warranty.


Correspondence: Marije ter Wal - m.j.terwal@bham.ac.uk

