function compile_evrything_on_the_cluster_mcc_VLP

cd ..;

global_path = pwd;

cd([global_path '/pipeline_VLP']);
system('mcc -m -R -nojvm -R -singleCompThread -v pipeline_VLP.m') ;
cd .. ;

cd([global_path '/pipeline_VLP_generate_trajectories']);
system('mcc -m -R -nojvm -R -singleCompThread -v pipeline_VLP_generate_trajectories.m') ;
cd .. ;

cd([global_path '/reload_vmesh_files_VLP_no_parallel']);
system('mcc -m -R -nojvm -R -singleCompThread -v reload_vmesh_files_VLP_no_parallel.m') ;
cd .. ;

cd([global_path '/spread_vmesh_to_origin_files_neurons_no_parallel']);
system('mcc -m -R -nojvm -R -singleCompThread -v spread_vmesh_to_origin_files_neurons_no_parallel.m') ;
cd .. ;

cd([global_path '/collect_all_clusters_files_VLP']);
system('mcc -m -R -nojvm -R -singleCompThread -v collect_all_clusters_files_VLP.m') ;
cd .. ;

cd([global_path '/full_post_treatment']);
system('mcc -m -R -nojvm -R -singleCompThread -v full_post_treatment.m') ;
cd .. ;




end