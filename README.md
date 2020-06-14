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

#### Run Ambuild
Paste the following into a file called grow.py
```
#!/usr/bin/env python3
import os, sys
ambuildHome = "/opt/ambuild"
paramsDir =  'params'
sys.path.insert(0, ambuildHome)

# This imports the builder cell module - this is the only module that should be required
from ambuild import ab_cell

# Create Cell and seed it with the blocks
cellDim=[50,50,50]
mycell = ab_cell.Cell(cellDim, atomMargin=0.5, bondMargin=0.5, bondAngleMargin=15, paramsDir=paramsDir)

#import the two fragment files if you have 2 different building blocks
fragA = "blocks/ch4.car"

mycell.libraryAddFragment( filename=fragA, fragmentType='A')
mycell.addBondType('A:a-A:a')

mycell.seed(1, fragmentType='A', center=True)
for i in range(1):
    mycell.growBlocks(toGrow=10, cellEndGroups='A:a', libraryEndGroups='A:a', maxTries=500)
    mycell.optimiseGeometry()
    mycell.dump()
```

Try running:

```/opt/ambuild/misc/ambuild_docker.py ./grow.py```
