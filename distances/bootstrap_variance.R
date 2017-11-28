library(MASS)

B=40;
N=80
p=0.2
prop=0.1
for (b in 1:B){
  ext=paste('boot',b,sep='')
  t=test_smooth_RD_changes(N,p,prop,T=11,verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",
                           path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_graph=ext)
}




N=80
spar=c(0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6)
prop=0.2
for (b in spar){
  ext=paste('spar_',b,sep='')
  t=test_smooth_RD_changes(N,b,prop,T=11,verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",
                           path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_graph=ext)
}


B=15;
m=30
N=60
for (b in 1:B){
  t1=test_smooth_Realistic_changes(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=1,loc=FALSE,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext=paste("Power_",b,sep=''),very_verbose=F)
  t2=test_smooth_Realistic_changes(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=2,loc=FALSE,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext=paste("Island_",b,sep=''),very_verbose=F)
  t3=test_smooth_Realistic_changes(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=3,loc=FALSE,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext=paste("DotProd_",b,sep=''),very_verbose=F)
  t4=test_smooth_Realistic_changes(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=4,loc=FALSE,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext=paste("SBM_",b,sep=''),very_verbose=F)
}


t=test_change_point(N=80,p=0.4,prop=0.05,prop2=0.2,p2=0.4, T=21,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="rdm")
t1=test_change_point_realistic(N=80,m=60,m2=60,p=0.2,p2=0.7,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=1, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=TRUE,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext="Power")
t2=test_change_point_realistic(N=80,m=60,m2=60,p=0.2,p2=0.7,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=2, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=TRUE,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext="Island")
t3=test_change_point_realistic(N=80,m=60,m2=60,p=0.2,p2=0.7,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=3, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=TRUE,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext="Dotprod")
t4=test_change_point_realistic(N=80,m=60,m2=60,p=0.2,p2=0.7,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=4, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=TRUE,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext="SBM")
    

test_smooth_RD_changes<-function(N,p,prop,T=11,verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_graph=''){
  