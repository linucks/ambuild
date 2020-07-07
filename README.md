# Amorphous Builder (Ambuild)
Ambuild is a python program for creating polymeric molecular structures.
Please feel free to follow Ambuild on [twitter](https://twitter.com/Ambuild2).

The code is developed by [Abbie Trewin's](https://twitter.com/AbbieTrewin) group at the [University of Lancaster](https://www.lancaster.ac.uk/sci-tech/about-us/people/abbie-trewin).

## Installation
In order to run at all, Ambuild requires [numpy](https://numpy.org/), which is easily installed into any Python installation with a command such as:
```
pip install numpy
```

With numpy installed Ambuild can be used to create molecular structures, but cannot run any Molecular Dynamics or Optimisation steps. In order to do that, [HOOMD-Blue](http://glotzerlab.engin.umich.edu/hoomd-blue/) is required, and in order to run systems of a reasonable size, HOOMD-Blue will need to be running on GPUs. If you already have HOOMD-Blue installed you are ready to start using Ambuild, otherwise the instructions below detail how to install Ambuild and HOOMD-Blue.

### Ubuntu/Debian

#### 1. Install Docker

Instructions from: https://docs.docker.com/engine/install/ubuntu/

1. Firstly, remove any old versions with the command:
  ```
sudo apt-get remove docker docker-engine docker.io containerd runc
```
  Update the list of packages and install those required to install Docker with the following two commands:
  ```
sudo apt-get update
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
```

2. Add Docker’s official GPG key (so that the downloaded packages can be validated) with the command:
```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```
3. Add the Docker 'stable' repository to the list of available repositories, so that apt can download from it:
```
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
```
4. With the Docker respository added to the list, update the list of packages and then install docker:
```
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

#### 2. Create Docker Group
Non-root users cannot run Docker by default, so it usally needs to be run under sudo; however this means any files created are owned by root, which is not a good idea. To allow users to run docker/ambuild without having sudo access, create docker group and add any users to it who will be required to run docker/ambuild.

1. First create the group for all users of docker. This may have already been done with the docker installation command, so it may not be required, but it's not a problem to run this command again
```
sudo groupadd docker
```

2. Add the current logged in user (specified by the $USER environment variable) to this group. Any other users can be added by replacing $USER in the below commmand with the unix username.
```
sudo usermod -aG docker $USER
````
3. Activate the group (this just saves logging out/in again)
```
newgrp docker
```
4. With Docker installed, the docker group set up and the current user added, test if you can run docker without using sudo:
```
docker run hello-world
```
If this works, you have a working Docker installation!

#### 3. Install NVIDIA Drivers
In order for applications within the Docker container to take advantage of GPU acceleration, you will need to install the NVIDIA GPU drivers for your card - the drivers are the piece of software that allow different programmes to communicate with the GPU card. There are instructions for how to do this on the [NVIDIA webssite](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)

On Ubunutu, the easiest way to do this seems to be with the command:
```
sudo ubuntu-drivers autoinstall
```

#### 4. Install NVIDIA Docker Runtime
Instructions from: https://github.com/NVIDIA/nvidia-docker

1. Run the following commands to install the nvidia-container-toolkit:

  ```
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

2. Test nvidia-smi with the latest official CUDA image

  ```docker run --gpus all nvidia/cuda:10.0-base nvidia-smi```

There is currently a bug with the nvidia docker container runtime as detailed here: https://github.com/docker/compose/issues/6691

To work around the bug carry out the following additional steps are required:

3. Install the nvidia-container-runtime
```
sudo apt install nvidia-container-runtime
```
4. Create a file called /etc/docker/daemon.json with the following content
```
{
    "runtimes": {
        "nvidia": {
            "path": "/usr/bin/nvidia-container-runtime",
            "runtimeArgs": []
        }
    }
}
```
  If you're unsure how to do this, cut and paste the following command into the terminal:
```
sudo tee -a /etc/docker/daemon.json << EOF
{
    "runtimes": {
        "nvidia": {
            "path": "/usr/bin/nvidia-container-runtime",
            "runtimeArgs": []
        }
    }
}
EOF
```
5. Restart docker:
```
sudo systemctl restart docker
```

#### 5. Disable secondary video GPU card
If you have more than one GPU card (e.g. you have a card specifically for running jobs), then you may need to disable your video GPU card for running jobs so that any GPU jobs are placed on the specialised card rather than using the video card. This will not disable the video card for viewing your screen - it will just prevent it being used to run computational simulation jobs.

1. Find ID of card to disable (this will print the UUID string, that you can then use in the command below).
```
nvidia-smi -L
```

2. Disable video GPU by setting mode to 2/PROHIBITED for the GPU card who's UUID we identified with the command above.
```
sudo nvidia-smi -c 2 -i GPU-4030396e-e7b4-aa4d-e035-22758536dba5
```

#### 6. Get Ambuild

1. Firstly, install git using:
```
sudo apt-get install git
```

2. Checkout ambuild from the github repository.
```
git clone https://github.com/linucks/ambuild.git
```

#### 7. Run Ambuild with Docker
To run Ambuild with docker a command like that below should be used. The key thing to understand is that the Docker container cannot _see_ the local computer filesystem - it can only access the directory structure within the container. In order to access files on the local computer, any directories will need to be mounted into the container using ```--volume``` arguments, and then the _internal_ container path used within any scripts.

> **NB: the backslash at the end of each line is a continuation character, so the whole block of text is actually a single command and could be typed as a single line.**

```
docker run --rm \
--runtime=nvidia \
--volume "$PWD":"$PWD" \
--workdir "$PWD" \
--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild \
--volume /home/abbie/Dropbox/Ambuild:/home/abbie/Dropbox/Ambuild \
glotzerlab/software \
python3 ambuild_script.py
```
Each line is explained below.

1. ```docker run --rm ``` run docker and remove the container on exit.

2. ```--runtime=nvidia``` use the Nvidia environment to take advantage of the GPU acceleration.
3. ```--volume "$PWD":"$PWD"``` make the current working directory from where this command is run (specified using the variable ```"$PWD"```) available inside the container.
4. ```--workdir $PWD``` make the working directory inside the container the full path to the working directory on the local machine. This means means that the current working directory where the script is run, will be used as the working directory for running Ambuild.
5. ```--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild``` make the directory ```/opt/ambuild/ambuild``` on the local filesystem available as ```/usr/lib/python3/dist-packages/ambuild``` within the container. This makes it possible for the python3 executable within the container to find the ambuild code, so that ```import ambuild``` within the ambuild_script.py works.
6. ```--volume /home/abbie/Dropbox/Ambuild:/home/abbie/Dropbox/Ambuild``` this is the only optional parameter. It makes the directory ```/home/abbie/Dropbox/Ambuild/params``` available within the container. Any additional paths to directories outside the current working directory that are used within the Ambuild script will need to be added in this way.
7. ```glotzerlab/software``` use the docker image from [glotzerlab/software](https://hub.docker.com/r/glotzerlab/software/). This downloads the file from the Docker repository (it's very large - several Gb - so the download can take some time, although it's only done once), and uses this to create the container.
8. ```python3 ambuild_script.py``` run the ```ambuild_script.py``` script, containing the Ambuild commands in the current directory with the python3 executable in the container.

##### 7.1 Run Ambuild with Docker using the run_ambuild_docker.sh script
The file [run_ambuild_docker.sh](https://github.com/linucks/ambuild/blob/master/misc/run_ambuild_docker.sh) that is distributed with Ambuild in the ```misc``` directory faciliates running Ambuild with an installed Docker installation. It creates the command to run the Docker container, and accepts additional ```--volume``` arguments, as well as the path to the ambuild script to run. An example of using it is below, where Ambuild has been downloaded and unpacked into the ```/opt/``` directory, and the files in the ```/home/abbie/Dropbox/Ambuild_Files``` directory need to be accessed within the ```cell_size_test.py``` Ambuild script:
```
/opt/ambuild/misc/run_ambuild_docker.sh \
--volume /home/abbie/Dropbox/Ambuild_Files:/home/abbie/Dropbox/Ambuild_Files \
cell_size_test.py
```
## Installation of optional dependencies
### Poreblazer
To install [poreblazer](https://github.com/richardjgowers/poreblazer) for use by Ambuild, the following steps are required.

1. Checkout or download poreblazer from github:
```git clone https://github.com/richardjgowers/poreblazer.git```

2. Install the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. On Ubuntu/Debian, this should just be a case of running:
```sudo apt-get install gfortran```

3. Compile the poreblazer executable. This is done in the ```src``` directory of the poreblazer directory, so cd into this directory and then run the command:```make``` This should create the ```poreblazer.exe``` executable in this directory.
