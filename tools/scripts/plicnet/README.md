
Conda
---------
1. Install [Anaconda](https://www.anaconda.com/download/) on your machine (make sure to select Python 3)

2. Follow the prompts in the terminal (you may need to run a conda shell script to fully define environment variables)

3. Make sure the installation went smoothly:
```bash
conda --version
```
If you get an error or don't see the version output, it's possible that your current bash doesn't know where conda is. This can be fixed by adding the anaconda or miniconda bin location to your systems PATH variable or to your *.bash_profile* or *.bashrc* 
```bash
# add path to conda
export PATH="Your_path_to_anaconda/anaconda3/bin:$PATH"
```
You may need to `source .bash_profile` or just restart your bash.

Creating a Virtual Envrionment
---------
1. Use Conda to create a virtual envrioment.
```bash
conda create -n myenv python=3.10
```

**NOTE: Python version >=3.7 but <= 3.11 is required for Pytorch**


2. If `myenv` is activated, you'll see `(myenv)` to the left of the terminal prompt. If not, you can activate it by
```bash
conda activate myenv
```

3. Onece `myenv` is activated, install the required python packages
```bash
pip install numpy
pip install torch
```

Running the Script
---------
1. Provide a Pytorch model file (.pt file) in the same directory as the generate_plicnet.py script. The script expects this file to be named "model.pt" by default. 

2. Run the script, and the resulting plicnet.f90 will be generated in the same directory. Copy plicnet.f90 to NGA2/src/two_phase