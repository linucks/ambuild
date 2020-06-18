# Amorphous Builder (Ambuild)
Ambuild is a python program for creating polymeric molecular structures.

## Installation
### Ubuntu/Debian

* install git
```sudo apt-get install git```

#### Install Docker
Instructions from: https://docs.docker.com/engine/install/ubuntu/
```
# Remove any old versions:
sudo apt-get remove docker docker-engine docker.io containerd runc
# Set up docker repository:
sudo apt-get update
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
# Add Docker’s official GPG key:
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
# Set up stable repository
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
# Install Docker
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```
To allow users to run docker/ambuild without having sudo access, create docker group and add any users to it who will be required to run docker/ambuild.

```
# NB: may not be required
sudo groupadd docker
# Add user abbie
sudo usermod -aG docker
# Active group (to save logging out/in)
newgrp docker
# Test can run hello world
docker run hello-world
```

#### Install NVIDIA Docker Runtime
Instructions: https://github.com/NVIDIA/nvidia-docker

```
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

Test nvidia-smi with the latest official CUDA image

```docker run --gpus all nvidia/cuda:10.0-base nvidia-smi```

Workaround for: https://github.com/docker/compose/issues/6691
```
# Seem to need this too:
sudo apt install nvidia-container-runtime
# Copy into /etc/docker/daemon.json
{
    "runtimes": {
        "nvidia": {
            "path": "/usr/bin/nvidia-container-runtime",
            "runtimeArgs": []
        }
    }
}
# Then:
sudo systemctl restart docker
```
**Disable video GPU.**
```
# Find ID of card to disable
nvidia-smi -L
# Disable video GPU by setting mode to 2/PROHIBITED
nvidia-smi -c 2 -i GPU-4030396e-e7b4-aa4d-e035-22758536dba5
```

### Checkout Ambuild
```
cd /opt
git clone https://github.com/linucks/ambuild.git
```

#### Run Ambuild with Docker
To run Ambuild with docker a command like that below should be used.

> **NB: the backslash at the end of each line is a continuation character, so the whole block of text is actually a single command and could be typed as a single line.**

The key thing to understand is that the Docker container cannot _see_ the local computer filesystem - it can only access the directory structure within the container. In order to access files on the local computer, any directories will need to be mounted into the container using ```--volume``` arguments, and then the _internal_ container path used within any scripts.

```
docker run -it --rm \
--runtime=nvidia \
--volume "$PWD":/home/glotzerlab \
--workdir /home/glotzerlab \
--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild \
--volume /home/abbie/Dropbox/Ambuild/params:/home/glotzerlab/params \
glotzerlab/software \
python3 ambuild_script.py
```
Each line is explained below.

1. ```docker run -it --rm ``` run docker in interactive mode and remove the container on exit.

2. ```--runtime=nvidia``` use the Nvidia environment to take advantage of the GPU acceleration.
3. ```--volume "$PWD":/home/glotzerlab``` Make the current working directory from where this command is run available inside the container as ```/home/glotzerlab```
4. ```--workdir /home/glotzerlab``` make the working directory inside the container ```/home/glotzerlab``` - this means that the current working directory where the script is run (specified using the variable ```"$PWD"```) - will be used as the working directory for running Ambuild.
5. ```--volume /opt/ambuild/ambuild:/usr/lib/python3/dist-packages/ambuild```
   >Make the directory ```/opt/ambuild/ambuild``` on the local filesystem available as ```/usr/lib/python3/dist-packages/ambuild``` within the container. This makes it possible for the python3 executable within the container to find the ambuild code, so that ```import ambuild``` within the ambuild_script.py works.
6. ```--volume /home/abbie/Dropbox/Ambuild/params:/home/glotzerlab/params```
   >This is the only optional parameter - it makes the directory ```/home/abbie/Dropbox/Ambuild/params``` available within the container as ```/home/glotzerlab/params```. Any additional paths to directories outside the current working directory path will need to be added in this way and the ***internal*** container path used within the Ambuild script i.e. the path to the params directory within the Ambuild script would be ```/home/glotzerlab/params``` **NOT** ```/home/abbie/Dropbox/Ambuild/params```.
7. ```glotzerlab/software```
   >Use the docker image from [glotzerlab/software](https://hub.docker.com/r/glotzerlab/software/). This downloads the file from the docker repository (it's very large - several Gb - so the download can take some time, although it's only done once), and uses this to create the container.
8. ```python3 ambuild_script.py``` Run the ```ambuild_script.py``` script, containing the Ambuild commands in the current directory with the python3 executable in the container.
