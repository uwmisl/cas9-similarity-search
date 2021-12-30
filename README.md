# cas9-similarity-search
DNA similarity search using Cas9

## A visual overview
You are currently in the `Cas9 Similarity Search` repository. 
The feature vectors we start with come from the preceding Open Images repository: https://github.com/uwmisl/primo-openimages 

![Image here](https://github.com/uwmisl/cas9-similarity-search/blob/main/documentation/similaritysearcharchitecture.png)


# Connecting to jupyterlab in a docker

## SSH Port Forwarding
Pick a port, P: anything unused; a large number.
When connecting with ssh: `ssh -L localhost:P:localhost:P brinstar.cs.washington.edu`. This forwards port P on your local workstation to port P on brinstar.

## Run docker
Launch the docker using the script in the cas9-similarity-search repository:
`sudo ./docker.sh -d <path to primo-openimages dir> -p P`, where P is the port you forwarded when connecting via SSH. The path is `/home/kstwrt/HDD/kendall/datasets/primo-openimages/`
`sudo ./docker.sh -d /home/kstwrt/HDD/kendall/datasets/primo-openimages/ -p P`

## Launch browser
Once you launch jupter notebook in the docker, you can connect in your browser. Copy the `/?token=<secret token>` bit from the jupyter output, and point your browser to:
`http://localhost:P/?token=<secret token>`.

## Launch shell in docker
`sudo docker exec -it <docker id> /bin/bash`
You can find your `docker id` by doing `sudo docker ps`

##  In docker, find Jupyter token to connect to
`jupyter notebook list`

# MISC Helpful Commands
For monitoring GPU:
`watch -n 1 nvidia-smi`