
[![Binder](https://cloud.githubusercontent.com/assets/15954923/21268015/28638b76-c3ad-11e6-89d5-6096e56f924a.png)](https://mybinder.org/v2/gh/naoufal51/Orion-Augmented/master?filepath=notebook.ipynb)
# Orion Augmented 
#### Work In Progress
In order to run the notebook showcasing the python version (partial) of the code used for the experiments, please use mybinder or clone the repository and use jupyter on you local machine for launching the notebook.

Please click on the orion logo to launch the notebook or here : [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/naoufal51/Orion-Augmented/master?filepath=notebook.ipynb)

We also provide code for CSI measurement using the experiment orchestration tool nepi-ng. It allows configuring the R2lab wireless nodes, run a transmission scenario between a transmitter and a receiver, and collect the CSI logs for angle of arrival (AoA) estimation. This code is tunable for more sophisticated transmission scenarios, for instance the usage of spatial multiplexing. https://github.com/parmentelat/r2lab-demos. We also propose an explanatory video for running an AoA experiment scenario https://r2lab.inria.fr/tuto-900-youtube.md#AOA.

# Orion Augmented: Orientation Estimation with Commodity Wi-Fi

With MIMO, Wi-Fi led the way to the adoption of antenna array signal processing techniques for fine-grained localization using commodity hardware. MIMO techniques, previously exclusive to specific domains of applications for instance radar systems, allow to consider estimating the orientation in space of a device equipped with COTS Wi-Fi chip. In fact, the availability channel state information measurement in recent Wi-Fi chips, makes it a candidate readily available for experimentation. Accordingly, we propose the ORION system to estimate the orientation (heading and yaw) of a MIMO Wi-Fi equipped object, relying on a joint estimation of the angle of arrival and the angle of departure.

<img src="/img/v0_lab.jpg" alt="Drawing" height="400px">
