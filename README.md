# RTW3D
Matlab codes for regularization of 3D teleseismic wavefield

Version: 2.0

Date: 06/30/2020

By

   Dr. Jinhai Zhang, zjh@mail.iggcas.ac.cn
   
   Dr. Chengliang Xie, clx@cugb.edu.cn

Acknowledgments:

   RBF function is from Matlab Central by Alex Chirokov (alex.chirokov@gmail.com). We also thank  C.D. Saragiotis (csar@auth.gr), Y.M. Altman (altmany@gmail.com) and F.J. Simons (fjsimons@alum.mit.edu) for their contributions on codes associated with the sac format files. 

Free for personal and non-commercial use.

We could not guarantee the codes work well all the time.

Thanks for feedback of any bug or suggestion.

===========================================================================

Main script: RTW3D.m

    Input:
    
        *.cfg: System parameters
        
            '#' defines comments;
            
            '>' defines specified information
            
        *.lst: file list of seismic data file in sac format, i.e. Observed stations
             
    Output:
    
         Output seismic stations in sac format, i.e. Regularized stations
             
===========================================================================

NOTES:

(1) Revision might be necessary based on different operating systems. The present version was conducted on linux.

(2) Several RBF functions are available for present version: 
    "Gaussian, Multiquadrics, Linear, Cubic and Thinplate"
    One could choose or redefine appropriate function. See comments in code.

(3) We assume that the waveforms have been aligned according to their P-wave arrival times.
    The input data could be pre-evaluated with function "preEval".

(4) For general cases, we set the smooth parameter as 0 for RBF method, and this could be user defined in RTW3D.m by searching 'smooth='

(5) For 3D RBF interpolation, the default setting is 'all in' (i.e. the entire dataset is loaded at one time). However, for general cases, this might be limited by available computer memory. Calculation has to be divided into 'local' 3D interpolations. The window length for each 'local dataset' could be defined in RTW3D.m by searching 'winLen='

===========================================================================

References:

[1] Buhmann, M.D. (2003), Radial Basis Functions: Theory and Implementations, Cambridge University Press.

[2] Fasshauer, G.E. (2007), Meshfree Approximation Methods with Matlab, World Scientific Publishers.

[3] Hu, S., X. Jiang, L. Zhu, and H. Yao (2019), Wavefield Reconstruction of Teleseismic Receiver Function with the Stretching‐and‐Squeezing Interpolation Method, Seismological Research Letters, 90(2A), 716-726, doi:10.1785/0220180197.

[4] Song, P., X. Zhang, Y. Liu, and J. Teng (2017), Moho imaging based on receiver function analysis with teleseismic wavefield reconstruction: Application to South China, Tectonophysics, 718, 118-131, doi:10.1016/j.tecto.2017.05.031.

[5] Zhang, J., and T. Zheng (2015), Receiver Function Imaging with Reconstructed Wavefields from Sparsely Scattered Stations, Seismological Research Letters, 86(1), 165-172, doi:10.1785/0220140028.

[6] Xie,C., Y. Fang, and J. Zhang (2021), Regularizing the 3D teleseismic wavefield for receiver function imaging using a radial basis function. Geophysical Journal International (In revision).

[7] Chirokov, A.(2006), Scattered Data Interpolation and Approximation using Radial Base Functions, edited, MATLAB Central File Exchange.
