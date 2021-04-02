# CANDO
## Computer Aided Nanostructure Design and Optimization

Cando is an open-source, computational chemistry environment for
building and designing functional molecules and materials.

<div style="text-align: center"><img src="documentation/CANDO_map.png" width=600 height=350 align=center/></div>

# Building Cando 
## Docker Image
The fastest method to start using CANDO is to download the [CANDO Docker image](https://hub.docker.com/r/thirdlaw/cando) from DockerHub or if you already have Docker installed:

```docker pull thirdlaw/cando```

Once downloaded, simply run:

```docker run -it -p 8888:8888 thirdlaw/cando```

and open a web-browser (Firefox is reccomended) at 0.0.0.0:8888

## Ubuntu 18.X
A CANDO binary can be downloaded from [ThirdLaw Technologies LLC](http://www.thirdlaw.tech/downloads/cando/binaries/ubuntu-18/).

To build CANDO on your native system, please see the state of dependencies for various installations on the [Wiki](https://github.com/cando-developers/cando/wiki).

# CANDO + Jupyter Demonstrations
Several JupyterLab notebooks that demonstrate how CANDO can be used in the realm of computational chemustry can be found in [CANDO-demos](https://github.com/cando-developers/cando-demos).

To run JupyterLab and start a demo:

1. run ```/opt/clasp/bin/cando-jupyter```

2. copy and paste the provided link into a web browser (Firefox is recommended).

# Video Tutorials
For up to date video demos and tutorials please [subscribe](https://www.youtube.com/channel/UC4xYBaHwB2kDCXaRALXdh7w) and like our videos on YouTube.

[Thermodynamic integration using a pre-build ligand library](https://www.youtube.com/watch?v=Md8rrrg7Kvg)

[Sketching molecules for thermodynamic integration calculations](https://www.youtube.com/watch?v=CSa2jzzwves)

[Running CANDO calculation jobs inside a Docker image](https://www.youtube.com/watch?v=3bMJB8dV2WU)

