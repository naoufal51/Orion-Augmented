![ScreenShot](https://cloud.githubusercontent.com/assets/15954923/21268015/28638b76-c3ad-11e6-89d5-6096e56f924a.png)
# Orion Augmented
In order to run the notebook showcasing the python version (partial) of the code used for the experiments, please use mybinder or clone the repository and use jupyter on you local machine for launching the notebook present in the code folder. If you choose to use mybinder, you should copy and paste the repository link into the appropriate field present in home page of the website: http://mybinder.org.

We also provide code for CSI measurement using the experiment orchestration tool nepi-ng. It allows configuring the R2lab wireless nodes, run a transmission scenario between a transmitter and a receiver, and collect the CSI logs for angle of arrival (AoA) estimation. This code is tunable for more sophisticated transmission scenarios, for instance the usage of spatial multiplexing. https://github.com/parmentelat/r2lab/tree/public/demos. We also propose an explanatory video for running an AoA experiment scenario https://r2lab.inria.fr/tuto-900-youtube.md#AOA.

# Orion Augmented Orientation Estimation Using Commodity Wi-Fi

With MIMO, Wi-Fi led the way to the adoption of antenna array signal processing techniques for fine-grained localization using commodity hardware. MIMO techniques, previously exclusive to specific domains of applications for instance radar systems, allow to consider estimating the orientation in space of a device equipped with COTS Wi-Fi chip. In fact, the availability channel state information measurement in recent Wi-Fi chips, makes it a candidate readily available for experimentation. Accordingly, we propose the ORION system to estimate the orientation (heading and yaw) of a MIMO Wi-Fi equipped object, relying on a joint estimation of the angle of arrival and the angle of departure.

<img src="/code/img/exp_setup2.jpg" alt="Drawing" height="400px">

