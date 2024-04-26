# Differentiable Voronoi Diagrams for Simulation of Cell-Based Mechanical Systems

Logan Numerow, [Yue Li](https://liyuesolo.github.io/), [Stelian Coros](https://crl.ethz.ch/people/coros/index.html), [Bernhard Thomaszewski](https://n.ethz.ch/~bthomasz/)

ACM Transactions on Graphics (Proc. SIGGRAPH 2024)

## Compile and Run 

For now we only provide support on Linux machines. The given instructions are for building with VSCode.

### Install NVIDIA Docker

> $ curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
    
> $ sed -i -e '/experimental/ s/^#//g' /etc/apt/sources.list.d/nvidia-container-toolkit.list

> $ sudo apt-get update \
    && sudo apt-get install -y nvidia-container-toolkit

> $ nvidia-ctk cdi generate --output=/etc/cdi/nvidia.yaml

> $ sudo nvidia-ctk runtime configure --runtime=docker

### Enable Display

Enable Docker to connect to host display to spawn GUI windows. Run from command line on host machine, repeat in case of display error.
> $ xhost +

### Run Docker in VSCode

Open the repository folder in VSCode. Install Docker and Dev Containers extensions in VSCode.

In VSCode, type `control + p`, then type `>Reopen in Container` (with the '>'). This option will show up in the >< tab in the bottom left corner of vscode.

The docker build will take a while (30 minutes).

That's it, you're ready to run the app!

## License
See the [LICENSE](./LICENSE) file for license rights and limitations.

## Contact
Please reach out to ([Logan Numerow](logan-numerow@hotmail.com)) for any issues or questions!

If this code contributes to your research, please consider citing our work:
```
@article{numerow2024voronoi,
  title={Differentiable Voronoi Diagrams for Simulation of Cell-Based Mechanical Systems},
  author={Numerow, Logan and Li, Yue and Coros, Stelian and Thomaszewski, Bernhard},
  journal={ACM Transactions on Graphics (TOG)},
  volume={43},
  number={4},
  year={2024},
  publisher={ACM New York, NY, USA}
}
```
