# CANDO
## Computer Aided Nanostructure Design and Optimization

Cando is an open-source, computational chemistry environment for
building and designing functional molecules and materials.

<div style="text-align: center"><img src="documentation/CANDO_map.png" width=600 height=350 align=center/></div>

## Installation

### Install Cando using package managers:

1. apt-get for Ubuntu 22.04, Debian Bookworm (in testing) and Windows Subsystem for Linux (WSL2).
2. homebrew for MacOS.
3. pacman for Arch Linux.

These installations can be automatically updated from the Nighly Package Repository for each operating system to get bug fixes and the most up-to-date version of Cando. We are fixing (and creating :-)) bugs all the time.

If you want to use the jupyter-lab kernel implemented in Cando then follow the entire instructions below for each operating system.  If you only want to install Cando - then skip the steps that mention jupyter-lab.

#### Installing from the Nightly Linux Package Repository.

Binary packages for Cando for Linux are built nightly starting about 2:00 AM EST and take about 4 hours to complete. Currently packages are built for Ubuntu 22.04, Debian Bookworm (in testing), and Arch. The Ubuntu 22.04 package works on WSL2 (Windows Subsystem for Linux). To add the package repository to your local system and install Cando execute one of the following commands in a shell. You may need to prefix these commands with `sudo` for root access.

1. Install jupyter-lab using miniconda.
  ```
3. Install Cando
  ```bash
  bash -c "$(curl -fsSL https://www.thirdlaw.tech/pkg/cando.sh)"
  ```
  
  If these scripts fail to detect your Operating System then you can install manually via one of the following.

- Debian Bookworm
    ```bash
    curl -fsSL https://www.thirdlaw.tech/pkg/key.gpg | gpg --yes --dearmor --output /etc/apt/trusted.gpg.d/thirdlaw.gpg
    add-apt-repository "deb https://www.thirdlaw.tech/pkg/bookworm/ ./"
    apt-get update
    apt-get install cando
    ```
- Ubunty Jammy
    ```bash
    curl -fsSL https://www.thirdlaw.tech/pkg/key.gpg | gpg --yes --dearmor --output /etc/apt/trusted.gpg.d/thirdlaw.gpg
    add-apt-repository "deb https://www.thirdlaw.tech/pkg/jammy/ ./"
    apt-get update
    apt-get install cando
    ```
- Arch
    ```bash
    curl -fsSL https://www.thirdlaw.tech/pkg/key.gpg | pacman-key --add -
    pacman-key --lsign-key 009C1F20737D81230A0928CBC812E10EC39DF202
    echo -e "\n[thirdlaw]\nServer = https://www.thirdlaw.tech/pkg/arch\n" >>/etc/pacman.conf
    pacman -Sy --noconfirm cando-git
    ```

3. Enhance the jupyter-lab/Cando environment by installing jupyter-lab widgets and visualization tools.
  ``` You may need to install nodejs to install 


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

