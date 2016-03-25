========
Overview
========

Motivations
-----------

Imagine the challenge of driving your car from home to 
work with your eyes closed. This is similar to how scientists 
feel today in conducting imaging experiments at the 
`Advanced Photon Source (APS) <http://aps.anl.gov>`_. They must 
set up their experiments by intuition, hoping to collect the 
maximum information about the material under study; however, 
they do not see the whole picture of what happened until 
sometime after the experiment is done. As a result, they do 
not know if they have modified the specimen by beam damage, 
or if they spent considerable time scanning uninteresting 
areas while collecting insufficient detail from the crucial 
features or a critical time point in a dynamic specimen. 

This package will let the driver see *while driving*, by  
algorithms to obtain preliminary images from selectively sampled 
data and using this for intelligent experimental control to 
collect the right amount of data from the right regions. A 
typical experiment in this new scheme will contain many loops; 
that is, the experiment will continuously inform the computation, 
and the computation will concurrently design the experiment.

Beneficeries
------------

`Tomography <https://en.wikipedia.org/wiki/Tomography>`_ 
is the premier technique at synchrotrons for studying 
thick samples, and is used with a large variety of contrast modes 
such as phase contrast, fluorescence and diffraction contrast. 
It has a long history of applications from materials research to 
biological, environmental and life sciences, and is one of the 
most utilized techniques at APS and at other synchrotrons worldwide.

Features
--------
* Generic geometry definitions for scanning probes
* Reconstruction algorithms for streaming data
* Visualization of geometries and reconstructions
