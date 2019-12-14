# DeepRF_SLR

DeepRF_SLR is multiband spin-echo pulse design method based on the SLR algorithm [1] and deep reinforcement learning. For more information, refer to the [submitted paper]()

## Method

![method_figure](/Figures/method_figure.png)

## Results

## Prerequisites

Before running DeepRF_SLR, you need to install,
1. Ubuntu 18.04.2 LTS
2. MATLAB 2019a
3. CVX package (http://cvxr.com/cvx/)
4. Python 3.6.8 with following libraries:
   - Tensorflow 1.13.1
   - Numpy 1.16.4
   - Scipy 1.3.0
5. MATLAB engine for python (https://www.mathworks.com/help/matlab/matlab-engine-for-python.html)

## Installing

0. Install all prerequisites
1. Download the repository and unzip it
2. Add directory path of unziped folder (with subfolder) to MATLAB search path

## Running

To run DeepRF_SLR, for example, type this commnad in Ubuntu terminal (in the directory where you unzip).
> python3 DeepRF_SLR.py 3 6 6 --gpu 0 --lr 1e-4 --nn 8 --nl 256 --gamma 1.0 --iter 100000

You can see the explanation of each argument using following command:
> python3 DeepRF_SLR.py --help  

After running the script, 'DeepRF_SLR_refo_design.mat' will be generated.

To see the design result, 
> run 'DeepRF_SLR_result_generation.m' in MATLAB.  

You can see the pulse shapes and simulated slice profile.  
'DeepRF_SLR_design_result.mat' will be generated, which contains many datas including RF pulses.

For comparison, we uploaded MATLAB script for the Monte-Carlo algorithm [2].  
> run 'Monte_Carlo_design.m' in MATLAB.  

This will generate 'MC_design_result.mat' that contains design results.

## License

This project is licensed under the MIT License - see the LICENSE file for details

## Acknowledgement

For implementation of DeepRF_SLR, we adopted and modifed MATLAB codes from:
1. [Root-flipped pulse design functions from Sharma's work](http://www.vuiis.vanderbilt.edu/~grissowa/) [2]
2. [Basic RF design functions from John Pauly's lab in Stanford](http://rsl.stanford.edu/research/software.html)
3. [Amplitude modulated pulse design functions from Seada's work](https://github.com/mriphysics/AM_multiband/) [3]

Thank you for all sharing softwares!

## References

1. J. Pauly, P. Le Roux, D. Nishimura, and A. Macovski, “Parameter relations for the Shinnar-Le Roux selective excitation pulse design algorithm,” IEEE Trans. Med. Imag., vol. 10, no. 1, pp. 53-65, Mar. 1991.
2. A. Sharma, M. Lustig, and W. A. Grissom, “Root‐flipped multiband refocusing pulses,” Magn. Reson. Med., vol. 75, no. 1, pp. 227-237, Jan. 2016.
3. S. Abo Seada, A. N. Price, J. V. Hajnal, and S. J. Malik, “Optimized amplitude modulated multiband RF pulse design,” Magn. Reson. Med., vol. 78, no. 6, pp. 2185-2193, Dec. 2017.

## Contact

Dongmyung Shin, Ph.D. candidate, Seoul National University.  
shinsae11@gmail.com  
http://list.snu.ac.kr
