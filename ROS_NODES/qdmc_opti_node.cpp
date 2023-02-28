/* Edoardo Caldarelli, Institut de Robòtica i Informàtica Industrial (CSIC-UPC)

code based on work by Adrià Luque Acera */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "general_functions.h"
#include "model_functions.h"
#include "mpc_functions.h"
#include <ros/ros.h>
#include "cartesian_msgs/CartesianCommand.h"
#include "geometry_msgs/Pose.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/Quaternion.h"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <mpc_pkg/TwoPoses.h>
#include <mpc_pkg/TwoPoints.h>
#include <mpc_pkg/OptiDataQDMC.h>
#include <mpc_pkg/HorizonControls.h>
#include <std_srvs/Empty.h>


using namespace std;
using namespace casadi;


// GLOBAL VARIABLES
string datapath = "/home/robot/Desktop/ALuque/mpc_ros/src/mpc_node/data/";
int NTraj = 0;
int nCOM = 0;
int Hp = 0;
int Hc = 0;
double Ts = 0.00;
bool shutdown_flag = false;
Eigen::MatrixXd in_params;
Eigen::Matrix3d Rcloth;
Eigen::VectorXd u_bef(6);
// ----------------

// MAIN CONTROL PARAMETERS
// --------------------------
// Objective Function Weights
double W_Q = 1*1e5;
double W_R = 1*1e7;
// Constraint bounds
double ubound = 0.02;
double gbound = 0.0;
// Structure variants
bool opt_du = false;
bool opt_Qa = false;
// --------------------------


// Linear model structure
typedef struct linmdl {
  int row;
  int col;
  double mass;
  double grav;
  double dt;
  Eigen::Vector3d stiffness;
  Eigen::Vector3d damping;
  double z_sum;
  Eigen::MatrixXd nodeInitial;
  Eigen::MatrixXd mat_x;
  Eigen::MatrixXd mat_y;
  Eigen::MatrixXd mat_z;
  Eigen::VectorXd coord_ctrl;
  Eigen::VectorXd coord_lc;
} LinMdl;


// ROS SUBSCRIBERS - CALLBACK FUNCTIONS
void inidataReceived(const mpc_pkg::OptiDataQDMC& msg) {
	
	// Data is:
	// - Vectors of "in_params", Traj_[]_Hp_rot, u_rot and x_ini_COM_rot, updated
	// - Rcloth, updated to the input data so new u_rot can be changed back
	// - u_bef, to add to the incremental u_lin and get absolute positions
	
	in_params.setZero(in_params.rows(), in_params.cols());

  in_params.block(0,0, 1, 6 * Hp)      = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                              (((vector<double>) msg.f_tilde).data(), 6 * Hp).transpose();
  in_params.block(1,0, 1,6) 	 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
  															 (((vector<double>) msg.u_rot).data(), 6).transpose();

	in_params.block(2,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_left_x).data(),  Hp+1).transpose();
	in_params.block(3,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_right_x).data(), Hp+1).transpose();
	in_params.block(4,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_left_y).data(),  Hp+1).transpose();
	in_params.block(5,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_right_y).data(), Hp+1).transpose();
	in_params.block(6,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_left_z).data(),  Hp+1).transpose();
	in_params.block(7,0, 1,Hp+1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                                 (((vector<double>) msg.traj_right_z).data(), Hp+1).transpose();
	//in_params.block(8,0, 3,Hp) = ... //d_hat
	
  Rcloth.col(0) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                  (((vector<double>) msg.Rcloth_x).data(), 3);
  Rcloth.col(1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                  (((vector<double>) msg.Rcloth_y).data(), 3);
  Rcloth.col(2) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
                  (((vector<double>) msg.Rcloth_z).data(), 3);
  
  u_bef = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
          (((vector<double>) msg.u_bef).data(), 6);
}


// ROS SERVICES - CALLBACK FUNCTIONS
bool shutdownCallback(std_srvs::Empty::Request &req, std_srvs::Empty::Response &resp) {
	shutdown_flag = true;
	ROS_INFO_STREAM("Shutting down Optimizer");
	return true;

}




// ---------------------------
// -------- MAIN PROG --------
// ---------------------------

int main(int argc, char **argv) {

	// Initialize the ROS system and become a node.
  ros::init(argc, argv, "mpc_opti_node");
  ros::NodeHandle rosnh;
  
	// Define client & server objects to all services
	ros::ServiceServer srv_shutdown = rosnh.advertiseService("node_shutdown", &shutdownCallback);
	
  //ros::ServiceClient clt_foo = rosnh.serviceClient<node::Service>("service_name");
  //ros::service::waitForService("service_name");
  
  // Define Publishers
  ros::Publisher pub_usomhp = rosnh.advertise<mpc_pkg::HorizonControls>
                              ("mpc_controller/u_som_hp", 1000);
  
  // Define Subscribers
  ros::Subscriber sub_inidata = rosnh.subscribe("mpc_controller/opti_inidata",
                                                 1000, &inidataReceived);
                                                 
  
  // Get parameters from launch
  if(!rosnh.getParam("/qdmc_opti_node/datapath", datapath)) {
  	ROS_ERROR("Need to define the datapath (where reference csv files are) parameter.");
  }
  if(!rosnh.getParam("/qdmc_opti_node/NTraj", NTraj)) {
  	ROS_ERROR("Need to define the NTraj (reference trajectory number) parameter.");
  }
  if(!rosnh.getParam("/qdmc_opti_node/nCOM", nCOM)) {
  	ROS_ERROR("Need to define the nCOM (COM mesh side size) parameter.");
  }
  if(!rosnh.getParam("/qdmc_opti_node/Hp", Hp)) {
  	ROS_ERROR("Need to define the Hp (prediction horizon) parameter.");
  }
  if(!rosnh.getParam("/qdmc_opti_node/Hc", Hc)) {
    ROS_ERROR("Need to define the Hc (control horizon) parameter.");
  }
  if(!rosnh.getParam("/qdmc_opti_node/Ts", Ts)) {
  	ROS_ERROR("Need to define the Ts (sample time) parameter.");
  }
  
  
  
 	// --------------------------------
	// 1. INITIAL PARAMETERS DEFINITION
	// --------------------------------
	
	// Load trajectories to follow

  Eigen::MatrixXd phi_l_Traj = getCSVcontent(datapath+"ref_"+to_string(NTraj)+"L.csv");
  Eigen::MatrixXd phi_r_Traj = getCSVcontent(datapath+"ref_"+to_string(NTraj)+"R.csv");
	
	// Get implied cloth size and position
	Eigen::MatrixXd dphi_corners1 = phi_r_Traj.row(0) - phi_l_Traj.row(0);
	double lCloth = dphi_corners1.norm();
  Eigen::Vector3d cCloth;
  cCloth = (phi_r_Traj.row(0) + phi_l_Traj.row(0))/2;
  cCloth(2) += lCloth/2;
  double aCloth = atan2(dphi_corners1(1), dphi_corners1(0));
  
  
  // Initialization of COM variable
  LinMdl COM;
  
  // Load parameter table and get row for COM
  // LUT cols: MdlSz, Ts, theta[7]
  Eigen::MatrixXd thetaLUT = getCSVcontent(datapath+"ThetaMdl_LUT.csv");
  bool COMparams = false;
  for (int LUTi=0; LUTi<thetaLUT.rows(); LUTi++) {
		if (thetaLUT(LUTi, 1) != Ts) continue;
		
		Eigen::MatrixXd theta_i = thetaLUT.block(LUTi,2, 1,7);
		
		// Set model parameters
		if (thetaLUT(LUTi, 0) == nCOM) {
			if (COMparams) ROS_WARN_STREAM("Several rows match the same COM parameters!");
			
		  COM.stiffness << theta_i(0), theta_i(1), theta_i(2);
  		COM.damping   << theta_i(3), theta_i(4), theta_i(5);
  		COM.z_sum      = theta_i(6);
  		COMparams = true;
		}
		
  }
  // If no row matched the settings
  if (!COMparams) {
  	ROS_WARN_STREAM("No rows match the COM parameters! Setting theta to default.");
  	
  	nCOM = 4;
		COM.stiffness << -200.0, -15.0, -200.0;
  	COM.damping   << -4.0, -2.5, -4.0;
  	COM.z_sum      =  0.03;
	}
  Eigen::MatrixXd S_tilde = getCSVcontent(datapath+"S_tilde.csv");
  // Rest of COM Definition
  int nxC  = nCOM;
  int nyC  = nCOM;
  COM.row  = nxC;
  COM.col  = nyC;
  COM.mass = 0.1;
  COM.grav = 9.8;
  COM.dt   = Ts;
  int COMlength = nxC*nyC;

	// Important Coordinates (upper/lower corners in x,y,z)
  Eigen::VectorXd COM_nodes_ctrl(2);
	Eigen::VectorXd COM_coord_ctrl(6);
	COM_nodes_ctrl << nyC*(nxC-1), nyC*nxC-1;
	COM_coord_ctrl << COM_nodes_ctrl(0), COM_nodes_ctrl(1),
										COM_nodes_ctrl(0)+nxC*nyC, COM_nodes_ctrl(1)+nxC*nyC,
										COM_nodes_ctrl(0)+2*nxC*nyC, COM_nodes_ctrl(1)+2*nxC*nyC;
	COM.coord_ctrl = COM_coord_ctrl;
	Eigen::VectorXd COM_coord_lc(6);
	COM_coord_lc << 0, nyC-1, nxC*nyC, nxC*nyC+nyC-1, 2*nxC*nyC, 2*nxC*nyC+nyC-1;
	COM.coord_lc = COM_coord_lc;
	
  // Define initial position of the nodes (for ext_force)
  // Second half of the vector is velocities (initial = 0)
  Eigen::MatrixXd posCOM(COMlength,3);
  posCOM = create_lin_mesh(lCloth, nCOM, cCloth, aCloth);
  
  Eigen::VectorXd x_ini_COM(2*3*COMlength);
  x_ini_COM.setZero(2*3*COMlength);
  
  x_ini_COM.segment(0, COMlength) = posCOM.col(0);
  x_ini_COM.segment(1*COMlength, COMlength) = posCOM.col(1);
  x_ini_COM.segment(2*COMlength, COMlength) = posCOM.col(2);
  
  // Rotate initial COM positions to XZ plane
  Eigen::Matrix3d RCloth_ini;
  RCloth_ini << cos(aCloth), -sin(aCloth), 0,
  				      sin(aCloth),  cos(aCloth), 0,
  				               0,             0, 1;               
	Eigen::MatrixXd posCOM_XZ(COMlength,3);
	posCOM_XZ = (RCloth_ini.inverse() * posCOM.transpose()).transpose();
  

 	// Get the linear model for the COM
 	Eigen::MatrixXd A(6*COMlength,6*COMlength);
 	Eigen::MatrixXd B(6*COMlength,6);
 	Eigen::VectorXd ext_force(6*COMlength);
 	A.setZero(6*COMlength,6*COMlength);
 	B.setZero(6*COMlength,6);
 	ext_force.setZero(6*COMlength);
 	
	tie(COM, A, B, ext_force) = init_linear_model(COM, posCOM_XZ);
  
	// ----------------------------------
	// 2. OPTIMIZATION PROBLEM DEFINITION
	// ----------------------------------

	// Solver timeout: less than prediction time
	double timeout_s = Ts*Hp/4;
	
	// Declare model variables
	int n_states = 2*3*COMlength;
	SX u = SX::sym("u", 6, Hc);
	
	// Convert eigen matrices to Casadi matrices
    
  DM S_tilde_DM = DM::zeros(6 * Hp, 6 * Hc);
  memcpy(S_tilde_DM.ptr(), S_tilde.data(), sizeof(double) * 6 * Hp * 6 * Hc);
	
	// Initial parameters of the optimization problem
	SX P = SX::sym("P", 1 + 1 + 6 + 3, 6 * Hp);
  SX f_tilde = P(0, Slice()).T();

  SX u0 = P(1,Slice(0,6)).T();
  SX Rp_mat = P(Slice(2,2+6), Slice(1, Hp + 1));
  SX Rp = Rp_mat(0, Slice()).T();
  
  for(int i=1; i<6; ++i){
      Rp = vertcat(Rp, Rp_mat(i, Slice()).T());
  }
	//SX d_hat = P(Slice(8,11), Slice(0,Hp));
  /*

	SX all_u = horzcat(u0, u);
	
	//ROS_INFO_STREAM(endl<<u_SOM<<endl);

	SX delta_u_mat = all_u(Slice(), Slice(1,Hc+1)) - all_u(Slice(), Slice(0,Hc));
  SX delta_u = delta_u_mat(0, Slice()).T();
  for(int i = 1; i < 6; ++i){
      delta_u = vertcat(delta_u, delta_u_mat(i, Slice()).T());
  }*/


	// Optimization variables
	SX w = u(0, Slice()).T();
	for (int i=1; i<6; i++) {
		w = vertcat(w, u(i,Slice()).T());
	}
	vector<double> lbw (6*Hc);
	vector<double> ubw (6*Hc);
	for (int i=0; i<6*Hc; i++) {
		lbw[i] = -ubound;
		ubw[i] = +ubound;
	}
	
	// Other variables of the opt problem
	SX obj_fun = 0.0;
	SX g;
	vector<double> ubg;
	vector<double> lbg;
	
	// Weights (adaptive) calculation: direction from actual to desired position
	SX Qa;
  // Pad initial output value, sliced by get_adaptive_Q
  SX x0 = DM::zeros(n_states, 1);
  for (int i = 1; i < 6; ++i){
      x0(COM_coord_lc(i)) = f_tilde(i * Hp);
  }

  if (opt_Qa == true) {
    Qa = get_adaptive_Q(Rp, x0, COM_coord_lc);
  }
  else {
    Qa = 1;
  }
  SX error = f_tilde + SX::mtimes(S_tilde_DM, w) - Rp;
  obj_fun = W_Q * SX::mtimes(SX::mtimes(error.T(), Qa), error) + SX::mtimes(SX::mtimes(w.T(), W_R), w);
  SX x_ctrl = u0;
  for (int k=0; k<Hc; k++) {              
    // Constraint: Constant distance between upper corners
    x_ctrl = x_ctrl + u(Slice(), k);
    g = vertcat(g, pow(x_ctrl(1) - x_ctrl(0), 2) + pow(x_ctrl(3) - x_ctrl(2), 2) + pow(x_ctrl(5) - x_ctrl(4), 2) - pow(lCloth,2));
    lbg.insert(lbg.end(), 1, -gbound);
    ubg.insert(ubg.end(), 1, +gbound);
  }
 

	// Encapsulate in controller object
	SXDict nlp_prob = {{"f",obj_fun},{"x",w},{"g",g},{"p",P}};
	Dict nlp_opts=Dict();
	nlp_opts["ipopt.print_level"] = 0;
	nlp_opts["ipopt.max_cpu_time"] = timeout_s;
	nlp_opts["ipopt.sb"] = "yes";
	nlp_opts["print_time"] = 0;

	Function controller = nlpsol("ctrl_sol", "ipopt", nlp_prob, nlp_opts);
	
	// ----------------------------------
	// 3. EXECUTION OF THE OPTIMIZER LOOP
	// ----------------------------------
	
	// Initial controls
	Eigen::VectorXd u_SOM(6);
	
	// Auxiliary variables for base changes and storage
	vector<double> u1_rotv;
	vector<double> u2_rotv;
	vector<double> u3_rotv;
	vector<double> u4_rotv;
	vector<double> u5_rotv;
	vector<double> u6_rotv;
  Eigen::VectorXd u1_rot;
  Eigen::VectorXd u2_rot;
  Eigen::VectorXd u3_rot;
  Eigen::VectorXd u4_rot;
  Eigen::VectorXd u5_rot;
  Eigen::VectorXd u6_rot;
	
	Eigen::MatrixXd uHp_rot(Hc,6);
	Eigen::MatrixXd uHp_p1_rot2(3,Hc);
	Eigen::MatrixXd uHp_p2_rot2(3,Hc);
	Eigen::MatrixXd uHp_p1_lin2(3,Hc);
	Eigen::MatrixXd uHp_p2_lin2(3,Hc);
	Eigen::MatrixXd uHp_lin(Hc,6);
	Eigen::MatrixXd uHp_SOM(Hc,6);
	
	
	// Resize input parameters matrix for all states
	in_params.resize(1 + 1 + 6 + 3, 6*Hp);
	
	// Wait for initial optidata
	ROS_INFO_STREAM("Initialized Optimizer");
	boost::shared_ptr<mpc_pkg::OptiDataQDMC const> OptiData0;
	OptiData0 = ros::topic::waitForMessage<mpc_pkg::OptiDataQDMC>("/mpc_controller/opti_inidata");
  
  
	// START LOOP
	ros::Rate rate(1/Ts);
	int index = 0;
	while(rosnh.ok() && !shutdown_flag) {
	
		// Save initial iteration time
		ros::Time iterT0 = ros::Time::now();
	
		// Check subscriptions (in_params, Rcloth, u_bef)
		ros::spinOnce();
	
		// Initial seed of the optimization by using the desired output changing rate
		Eigen::MatrixXd dRef = in_params.block(2, 1, 6, Hc) - in_params.block(2, 0, 6, Hc);
		Eigen::Map<Eigen::VectorXd> args_x0(dRef.data(), dRef.size());
	
		// Transform variables for solver
		DM x0_dm = DM::zeros(6*Hc, 1);
		memcpy(x0_dm.ptr(), args_x0.data(), sizeof(double) * 6 * Hc);
	
		DM p_dm = DM::zeros(in_params.rows(), in_params.cols());
		memcpy(p_dm.ptr(), in_params.data(), sizeof(double)*in_params.rows()*in_params.cols());
	
		// Create the structure of parameters
		map<string, DM> arg, sol;
		arg["lbx"] = lbw;
		arg["ubx"] = ubw;
		arg["lbg"] = lbg;
		arg["ubg"] = ubg;
		arg["x0"]  = x0_dm;
		arg["p"]   = p_dm;
	
		// Find the solution
		sol = controller(arg);
		if(controller.stats()["return_status"]!="Solve_Succeeded"){
		  ROS_ERROR_STREAM(controller.stats()["return_status"] << std::endl);
		}
		//if(index==3) return 0;
		// Check how long it took
		ros::Duration optiDT = ros::Time::now() - iterT0;
		//int optiSteps = ceil(optiDT.toSec()/Ts);
		//ROS_INFO_STREAM("Opt.time: "<<1000*optiDT.toSec()<<" ms \t("<<optiSteps<<" steps)");
		if (optiDT.toSec() >= timeout_s) {
			int optiSteps = ceil(optiDT.toSec()/Ts);
		  ROS_WARN_STREAM("SOLVER TIMED OUT ("<<
		                  1000*optiDT.toSec() <<" ms / "<<
		                  optiSteps<<" steps)");
		  continue;
		}
		

		// Get control actions from the solution
		//  They are upper corner displacements (incremental pos)
		//  And they are in local cloth base (rot)
		DM wsol = sol["x"];
		//std:: cout << wsol << std::endl;
		//sleep(8);
		DM usol = DM::zeros(Hc,6);
		for (int i=0; i<6; i++) {
			usol(Slice(), i) = wsol(Slice(i*Hc, (i+1)*Hc));
		}
		//if (index%100==0){		
		// ROS_WARN_STREAM(usol);
		// sleep(10);
		//}
		
		// Process controls for the whole horizon
		u1_rotv = (vector<double>) usol(Slice(),0); //x1
		u2_rotv = (vector<double>) usol(Slice(),1); //x2
		u3_rotv = (vector<double>) usol(Slice(),2); //y1
		u4_rotv = (vector<double>) usol(Slice(),3); //y2
		u5_rotv = (vector<double>) usol(Slice(),4); //z1
		u6_rotv = (vector<double>) usol(Slice(),5); //z2
	  u1_rot = Eigen::VectorXd::Map(u1_rotv.data(), u1_rotv.size());
	  u2_rot = Eigen::VectorXd::Map(u2_rotv.data(), u2_rotv.size());
	  u3_rot = Eigen::VectorXd::Map(u3_rotv.data(), u3_rotv.size());
	  u4_rot = Eigen::VectorXd::Map(u4_rotv.data(), u4_rotv.size());
	  u5_rot = Eigen::VectorXd::Map(u5_rotv.data(), u5_rotv.size());
	  u6_rot = Eigen::VectorXd::Map(u6_rotv.data(), u6_rotv.size());
	  uHp_rot << u1_rot, u2_rot, u3_rot, u4_rot, u5_rot, u6_rot;
	  
	  uHp_p1_rot2.row(0) = u1_rot.transpose();
	  uHp_p1_rot2.row(1) = u3_rot.transpose();
	  uHp_p1_rot2.row(2) = u5_rot.transpose();
	  uHp_p2_rot2.row(0) = u2_rot.transpose();
	  uHp_p2_rot2.row(1) = u4_rot.transpose();
	  uHp_p2_rot2.row(2) = u6_rot.transpose();
	  
	  uHp_p1_lin2 = Rcloth * uHp_p1_rot2;
	  uHp_p2_lin2 = Rcloth * uHp_p2_rot2;
	  
	  uHp_lin.col(0) = uHp_p1_lin2.row(0).transpose();
	  uHp_lin.col(1) = uHp_p2_lin2.row(0).transpose();
	  uHp_lin.col(2) = uHp_p1_lin2.row(1).transpose();
	  uHp_lin.col(3) = uHp_p2_lin2.row(1).transpose();
	  uHp_lin.col(4) = uHp_p1_lin2.row(2).transpose();
	  uHp_lin.col(5) = uHp_p2_lin2.row(2).transpose();
	  
	  // u_SOM(n) = u_bef + Sum(u_lin(i))_(i=0 to n-1)
	  uHp_SOM.row(0) = uHp_lin.row(0) + u_bef.transpose();
	  for (int i=1; i<Hc; i++) {
	  	uHp_SOM.row(i) = uHp_lin.row(i) + uHp_SOM.row(i-1);
	  }
	  
		// Publish control actions
		mpc_pkg::HorizonControls uHp_SOM_pub;
		uHp_SOM_pub.u1Hp = vector<double>(uHp_SOM.col(0).data(),
								uHp_SOM.col(0).size()+uHp_SOM.col(0).data());
		uHp_SOM_pub.u2Hp = vector<double>(uHp_SOM.col(1).data(),
								uHp_SOM.col(1).size()+uHp_SOM.col(1).data());
		uHp_SOM_pub.u3Hp = vector<double>(uHp_SOM.col(2).data(),
								uHp_SOM.col(2).size()+uHp_SOM.col(2).data());
		uHp_SOM_pub.u4Hp = vector<double>(uHp_SOM.col(3).data(),
								uHp_SOM.col(3).size()+uHp_SOM.col(3).data());
		uHp_SOM_pub.u5Hp = vector<double>(uHp_SOM.col(4).data(),
								uHp_SOM.col(4).size()+uHp_SOM.col(4).data());
		uHp_SOM_pub.u6Hp = vector<double>(uHp_SOM.col(5).data(),
								uHp_SOM.col(5).size()+uHp_SOM.col(5).data());
	  
	  pub_usomhp.publish(uHp_SOM_pub);		
		// Debugging
	  //ROS_INFO_STREAM(endl<<u_SOM<<endl);
		
		// Execute at a fixed rate
		rate.sleep();
		++index;
	
	}
	// END LOOP
	
  
}































