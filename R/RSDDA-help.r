# I had to pull the following documentation from the help folder in the RSDDA package.
# I could not simply load the package with library(RSDDA) because I received the error message:
#	Error: package 'RSDDA' was built before R 2.10.0: please re-install it

RSDDA                 package:RSDDA                 R Documentation

_R_e_g_u_l_a_r_i_z_e_d _S_h_r_i_n_k_a_g_e-_b_a_s_e_d _D_i_a_g_o_n_a_l _D_i_s_c_r_i_m_i_n_a_n_t _A_n_a_l_y_s_i_s

_D_e_s_c_r_i_p_t_i_o_n:

     Regularized Shrinkage-based Diagonal Discriminant Analysis.

_U_s_a_g_e:

     RSDDA(infile, genenum, sampsize, topgenes, numofruns)

_A_r_g_u_m_e_n_t_s:

  infile: a data file with genes in rows and subjects in columns -
          class 1 followed by class 2.

 genenum: number of genes in the data file.

sampsize: sample size in class 1 or class 2 - balanced case.

topgenes: top genes chosen using the ratio of between-group to within
          group sums of square.

numofruns: number of runs.

_V_a_l_u_e:

mean misclassification error: Mean misclassification error of RSDDA
          over the number of runs.

_A_u_t_h_o_r(_s):

     Herbert Pang <pathwayRF@gmail.com>, Tiejun Tong
     <tiejun.tong@colorado.edu>, and Hongyu Zhao <hongyu.zhao@yale.edu>

_R_e_f_e_r_e_n_c_e_s:

_o "Regularized Shrinkage-based Diagonal Discriminant Analysis. 
     Biometrics (submitted) <URL:
     http://bioinformatics.med.yale.edu/rsdda/rsdda.htm>

_S_e_e _A_l_s_o:

     'MASS'.

_E_x_a_m_p_l_e_s:

     ## load example data
     data(sm8vs8)
     write.table(sm8vs8, file="sm8vs8.txt", sep="\t", row.names=FALSE, col.names=FALSE)
     ## run RSDDA with 10 runs
     ## using hypothetical data set with 2725 genes
     ## class 1 = class 2 = 8 samples
     ## select top 100 genes
     ## output mean misclassification rate
     RSDDA("sm8vs8.txt", genenum=2725, sampsize=8, topgenes=100, numofruns=10)

